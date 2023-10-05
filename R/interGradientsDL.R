#' @importFrom tensorflow %as%
NULL

################################################################################
############### Calculation of gradients for NN interpretation #################
################################################################################

#' Calculate gradients with respect to predicted cell types/loss for 
#' interpreting a trained deconvolution model
#'
#' Train a deep neural network model using training data from the
#' \code{\linkS4class{SpatialDDLS}} object. In addition, the trained model is
#' evaluated using test data, and prediction results are obtained to determine
#' its performance (see \code{?\link{calculateEvalMetrics}}).
#'
#' @param object \code{\linkS4class{SpatialDDLS}} object with a trained 
#'    deconvolution model (\code{trained.model} slot) and pure mixed 
#'    transcriptional profiles (\code{mixed.profiles} slot).
#' @param method Method to calculate gradients with respect to inputs. It can be
#'    \code{'class'} (gradients of predicted classes w.r.t. inputs), 
#'    \code{'loss'} (gradients of loss w.r.t. inputs) or \code{'both'}.
#' @param normalize Whether normalize data using logCPM (\code{TRUE} by 
#'   default). This parameter is only considered when the method used to 
#'   simulate the mixed transcriptional profiles (\code{simMixedProfiles} 
#'   function) was \code{"AddRawCount"}. Otherwise, data were already 
#'   normalized. This parameter should be set according to the transformation 
#'   used to train the model. 
#' @param scaling How to scale data. It can be:
#'   \code{"standarize"} (values are centered around the mean with a unit
#'   standard deviation), \code{"rescale"} (values are shifted and rescaled so
#'   that they end up ranging between 0 and 1, by default) or \code{"none"} (no
#'   scaling is performed). This parameter should be set according to the 
#'   transformation used to train the model. 
#' @param on.the.fly Not available.
#' @param agg.function Function used to build mixed transcriptional profiles. It
#'   may be: \itemize{ \item \code{"AddRawCount"} (by default): single-cell
#'   profiles (raw counts) are added up across cells. Then, log-CPMs are
#'   calculated. \item \code{"MeanCPM"}: single-cell profiles (raw counts) are
#'   transformed into CPMs and cross-cell averages are calculated. Then,
#'   \code{log2(CPM + 1)} is calculated. \item \code{"AddCPM"}: single-cell
#'   profiles (raw counts) are transformed into CPMs and are added up across
#'   cells. Then, log-CPMs are calculated.}
#' @param threads Number of threads used during simulation of mixed
#'   transcriptional profiles if \code{on.the.fly = TRUE} (1 by default).
#' @param verbose Show informative messages during the execution (\code{TRUE} by
#'   default).
#'
#' @return Object containing gradients in the \code{interpret.gradients} slot of
#'   the \code{DeconvDLModel} object (\code{trained.model} slot).
#'
#' @export
#'
#' @seealso \code{\link{plotTrainingHistory}} \code{\link{deconvSpatialDDLS}}
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(
#'     counts = matrix(
#'       rpois(30, lambda = 5), nrow = 15, ncol = 10,
#'       dimnames = list(paste0("Gene", seq(15)), paste0("RHC", seq(10)))
#'     )
#'   ),
#'   colData = data.frame(
#'     Cell_ID = paste0("RHC", seq(10)),
#'     Cell_Type = sample(x = paste0("CellType", seq(2)), size = 10,
#'                        replace = TRUE)
#'   ),
#'   rowData = data.frame(
#'     Gene_ID = paste0("Gene", seq(15))
#'   )
#' )
#' SDDLS <- createSpatialDDLSobject(
#'   sc.data = sce,
#'   sc.cell.ID.column = "Cell_ID",
#'   sc.gene.ID.column = "Gene_ID",
#' )
#' SDDLS <- genMixedCellProp(
#'   object = SDDLS,
#'   cell.ID.column = "Cell_ID",
#'   cell.type.column = "Cell_Type",
#'   num.sim.spots = 50,
#'   verbose = TRUE
#' )
#' SDDLS <- simMixedProfiles(SDDLS)
#' SDDLS <- trainDeconvModel(
#'   object = SDDLS,
#'   batch.size = 12,
#'   num.epochs = 5
#' )
#' ## calculating gradients
#' SDDLS <- interGradientsDL(SDDLS)
#' }
#'   
interGradientsDL <- function(
    object,
    method = "class",
    normalize = TRUE,
    scaling = "rescale",
    on.the.fly = FALSE,
    agg.function = "AddRawCount",
    threads = 1,
    verbose = TRUE
) {
  if (!is(object, "SpatialDDLS")) {
    stop("The provided object is not of SpatialDDLS class")
  } else if (is.null(trained.model(object))) {
    stop("The provided object does not have a trained model for evaluation")
  } else if (is.null(prob.cell.types(object)) ||
             !c("train", "test") %in% names(prob.cell.types(object))) {
    stop("The provided object does not contain actual cell proportions in ", 
         "'prob.cell.types' slot for test data")
  }
  if (!scaling %in% c("standarize", "rescale", "none")) {
    stop("'scaling' argument must be one of the following options: 'standarize', 'rescale' or 'none'")
  } else {
    if (scaling == "standarize") {
      scaling.fun <- base::scale
    } else if (scaling == "rescale") {
      scaling.fun <- rescale.function
    } else if (scaling == "none") {
      scaling.fun <- function(x) return(x)
    }
  }
  ## check if in prob cell types there are pure mixed transcriptional profiles
  num.pure <- rbind(
    prob.cell.types(object, type.data = "train") %>% prob.matrix(),
    prob.cell.types(object, type.data = "test") %>% prob.matrix()
  ) 
  num.pure <- colSums(num.pure == 100)
  if (any(num.pure == 0)) {
    stop(
      paste0(
        "Not all cell types have pure mixed transcriptional profiles, i.e. ", 
        "transcriptional profiles made of only one cell type. See ?genMixedCellProp", 
        " to generate them"
      )
    )
  }
  
  if (!on.the.fly) {
    ## checking if all data are provided
    if (is.null(mixed.profiles(object, type.data = "train")) | 
        is.null(mixed.profiles(object, type.data = "test"))) {
      stop("If on.the.fly == FALSE, mixed transcriptional profiles must be provided")
    } 
    ## checking number of pure mixed transcriptional profiles 
    num.pure <- rbind(
      as.matrix(mixed.profiles(object, type.data = "train")@colData),
      as.matrix(mixed.profiles(object, type.data = "test")@colData)
    ) 
    num.pure <- colSums(num.pure == 100)    
    if (any(num.pure == 0)) {
      stop(
        paste0(
          "Not all cell types have pure mixed transcriptional profiles, i.e. ", 
          "transcriptional profiles made of only one cell type. See ?genMixedCellProp", 
          " and ?simMixedProfiles to generate them"
        )
      )
    }
    ## get data and metadata 
    metadata.prop <- rbind(
      as.matrix(mixed.profiles(object, type.data = "train")@colData),
      as.matrix(mixed.profiles(object, type.data = "test")@colData)
    )
    metadata.prop <- metadata.prop[apply(X = metadata.prop == 100, MARGIN = 1, FUN = any), ]
    data <- cbind(
      assays(mixed.profiles(object, type.data = "train"))[["counts"]],
      assays(mixed.profiles(object, type.data = "test"))[["counts"]]
    )[, rownames(metadata.prop)]
    ## normalization with logCPM (if required)
    mixing.fun <- mixed.profiles(object, type.data = "train")@metadata[["mixing.fun"]]
    if (normalize & mixing.fun == "AddRawCount") {
      data <- log2(.cpmCalculate(x = data + 1))
    }
    ## standarization
    data <- scaling.fun(t(data))
  } else {
    stop("This functionality is not available yet. Try on.the.fly = FALSE")
    if (verbose) message("=== Generating mixed transcriptional profiles on the fly")
    ## just in case of on.the.fly = TRUE
    if (!agg.function %in% c("MeanCPM", "AddCPM", "AddRawCount")) {
      stop("'agg.function' must be one of the following options: 'MeanCPM', 'AddCPM', 'AddRawCount'")
    } else {
      if (agg.function == "MeanCPM") {
        .agg.fun <- aggregation.fun.mean.cpm
      } else if (agg.function == "AddCPM") {
        .agg.fun <- aggregation.fun.add.cpm
      } else if (agg.function == "AddRawCount") {
        .agg.fun <- aggregation.fun.add.raw.counts
      }
    }
    ## generate the matrices: not implemented yet
  }
  # checking if DNN model is json format or compiled
  if (is.list(trained.model(object)@model)) {
    model.comp <- .loadModelFromJSON(trained.model(object))
    trained.model(object) <- model.comp
  }
  if (method == "class") {
    gradients <- list(
      class = .calcGradientsClass(
        x.data = data, y.metadata = metadata.prop, model = trained.model(object)@model
      )
    )
  } else if (method == "loss") {
    gradients <- list(
      loss = .calcGradientsLoss(
        x.data = data, y.metadata = metadata.prop, model = trained.model(object)@model
      )
    )
  } else {
    gradients <- list(
      class = .calcGradientsClass(
        x.data = data, y.metadata = metadata.prop, model = trained.model(object)@model
      ),
      loss = .calcGradientsLoss(
        x.data = data, y.metadata = metadata.prop, model = trained.model(object)@model
      )
    )
  }
  object@trained.model@interpret.gradients <- gradients
  
  return(object)
}

## core function: this function is intended to receive the data and metadata 
## of those samples used to calculate gradients wrt class. previous functions
## will provide with the correct samples
.calcGradientsClass <- function(x.data, y.metadata, model) { 
  ## info
  n.samples <- nrow(x.data)
  samples.name <- rownames(x.data)
  features <- colnames(x.data)
  cell.types <- colnames(y.metadata)
  ## from matrix to tensorflow
  y.metadata <- y.metadata / 100
  x.data <- tensorflow::tf$Variable(x.data)  # Convert all samples at once
  y.metadata.tf <- tensorflow::tf$Variable(as.matrix(y.metadata))
  # Define the gradient tape
  list.gradients <- list()
  with(tensorflow::tf$GradientTape(persistent = TRUE) %as% tape, {
    # Forward pass through the model for all samples
    y.metadata.pred <- model(x.data)
    for (i in seq_along(cell.types)) {
      list.gradients[[cell.types[i]]] <- tape$gradient(
        y.metadata.pred[, i], x.data
      )
    }
  })
  gradients.matrix <- lapply(
    names(list.gradients), \(i) {
      gradients <- list.gradients[[i]]  %>% as.matrix()
      rownames(gradients) <- samples.name
      colnames(gradients) <- features
      return(gradients[y.metadata[, i] == 1, , drop = FALSE])
    }
  )  %>% do.call(rbind, .)
  
  return(gradients.matrix[rownames(y.metadata), , drop = FALSE])
}

## vanilla gradient is defined as the gradient of the loss function for the class
## we are interested in wrt the input variables
.calcGradientsLoss <- function(x.data, y.metadata, model) {
  ## info
  n.samples <- nrow(x.data)
  samples.name <- rownames(x.data)
  features <- colnames(x.data)
  cell.types <- colnames(y.metadata)
  ## from matrix to tensorflow
  y.metadata <- y.metadata / 100
  x.data <- tensorflow::tf$Variable(x.data)  # Convert all samples at once
  y.metadata.tf <- tensorflow::tf$Variable(as.matrix(y.metadata))
  # Define the gradient tape
  list.gradients <- list()
  with(tensorflow::tf$GradientTape(persistent = TRUE) %as% tape, {
    # Forward pass through the model for all samples
    y.metadata.pred <- model(x.data)
    loss <- keras::metric_kullback_leibler_divergence(y.metadata, y.metadata.pred) 
    gradients <- tape$gradient(loss, x.data)
  })
  gradients.matrix <- as.matrix(gradients)
  colnames(gradients.matrix) <- features
  rownames(gradients.matrix) <- samples.name
  
  return(gradients.matrix)
}


################################################################################
######################### Top gradients per cell type ##########################
################################################################################

#' Get top genes with largest/smallest gradients per cell type
#'
#' Get feature names of top features with largest/smallest gradients per cell type. 
#' These genes can be used as input to visualize how they are spatially expressed 
#' (\code{plotGeneSpatial} function) or to plot gradients a a heatmap (\code{plotGradHeatmap} function).
#'
#' @param object \code{\linkS4class{SpatialDDLS}} object with a trained 
#'    deconvolution model (\code{trained.model} slot) and pure mixed 
#'    transcriptional profiles (\code{mixed.profiles} slot).
#' @param method Method to calculate gradients with respect to inputs. It can be
#'    \code{'class'} (gradients of predicted classes w.r.t. inputs), 
#'    \code{'loss'} (gradients of loss w.r.t. inputs) or \code{'both'}. 
#'    \code{'class'} by default.
#' @param top.n Top n genes (positive and negative) taken per cell type. 
#'
#' @return Object containing gradients in the \code{interpret.gradients} slot of
#'   the \code{DeconvDLModel} object (\code{trained.model} slot).
#'
#' @export
#'
#' @seealso \code{\link{interGradientsDL}} \code{\link{trainDeconvModel}}
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(
#'     counts = matrix(
#'       rpois(30, lambda = 5), nrow = 15, ncol = 10,
#'       dimnames = list(paste0("Gene", seq(15)), paste0("RHC", seq(10)))
#'     )
#'   ),
#'   colData = data.frame(
#'     Cell_ID = paste0("RHC", seq(10)),
#'     Cell_Type = sample(x = paste0("CellType", seq(2)), size = 10,
#'                        replace = TRUE)
#'   ),
#'   rowData = data.frame(
#'     Gene_ID = paste0("Gene", seq(15))
#'   )
#' )
#' SDDLS <- createSpatialDDLSobject(
#'   sc.data = sce,
#'   sc.cell.ID.column = "Cell_ID",
#'   sc.gene.ID.column = "Gene_ID",
#' )
#' SDDLS <- genMixedCellProp(
#'   object = SDDLS,
#'   cell.ID.column = "Cell_ID",
#'   cell.type.column = "Cell_Type",
#'   num.sim.spots = 50,
#'   verbose = TRUE
#' )
#' SDDLS <- simMixedProfiles(SDDLS)
#' SDDLS <- trainDeconvModel(
#'   object = SDDLS,
#'   batch.size = 12,
#'   num.epochs = 5
#' )
#' ## calculating gradients
#' SDDLS <- interGradientsDL(SDDLS)
#' listGradients <- topGradientsCellType(SDDLS)
#' lapply(listGradients, head, n = 5)
#' }
#'   
topGradientsCellType <- function(object, method = "class", top.n.genes = 15) {
  if (!is(object, "SpatialDDLS")) {
    stop("The provided object is not of SpatialDDLS class")
  } else if (is.null(trained.model(object))) {
    stop("The provided object does not have a trained model for evaluation")
  }
  ## check that the number of top.n.genes is less than the number of genes
  if (length(trained.model(object)@features) < top.n.genes) {
    stop("top.n.genes argument is too large for the number of genes used to train the model. Set a lower number.")
  }
  ## get metadata
  metadata.prop <- rbind(
    as.matrix(mixed.profiles(object, type.data = "train")@colData),
    as.matrix(mixed.profiles(object, type.data = "test")@colData)
  )
  metadata.prop <- metadata.prop[apply(X = metadata.prop == 100, MARGIN = 1, FUN = any), ]
  if (method == "both") {
    ## check if both are present in the object
    if (!all(c("class", "loss") %in% names(trained.model(object)@interpret.gradients))) {
      stop("If method == 'both', 'class' and 'loss' gradients must be present in the SpatialDDLS object")
    }
    grads.top <- lapply(
      c("class", "loss"), \(met) {
        grad <- trained.model(object)@interpret.gradients[[met]]
        return(top.gradients(grad = grad, metadata = metadata.prop, n = top.n.genes))
      }
    ) %>% setNames(c("class", "loss"))
  } else if (method == "class" | method == "loss") {
    if (!any(method %in% names(trained.model(object)@interpret.gradients))) {
      stop("Chosen method is not present in the SpatialDDLS bject")
    }
    grad <- trained.model(object)@interpret.gradients[[method]]
    grads.top <- top.gradients(grad = grad, metadata = metadata.prop, n = top.n.genes)
  } else {
    stop("method parameter has to be one of the following options: 'both', 'class' or 'loss'")
  }

  return(grads.top)
}

top.gradients <- function(grad, metadata, n) {
  lapply(
    X = colnames(metadata), \(x) {
      spots <- rownames(metadata)[metadata[, x] == 100]
      mean.grads <- base::sort(x = colMeans(grad[spots, , drop = FALSE]), decreasing = TRUE)
      return(
        list(
          Positive = names(head(mean.grads, n = n)), 
          Negative = names(tail(mean.grads, n = n))
        )
      )
    }
  ) %>% setNames(colnames(metadata))
}


################################################################################
############################ Heatmap of gradients ##############################
################################################################################

#' Plot a heatmap of gradients
#'
#' Get feature names of top features with largest/smallest gradients per cell type. 
#' These genes can be used as input to visualize how they are spatially expressed 
#' (\code{plotGeneSpatial} function) or to plot gradients a a heatmap (\code{plotGradHeatmap} function).
#'
#' @param object \code{\linkS4class{SpatialDDLS}} object with a trained 
#'    deconvolution model (\code{trained.model} slot) and pure mixed 
#'    transcriptional profiles (\code{mixed.profiles} slot).
#' @param method Method to calculate gradients with respect to inputs. It can be
#'    \code{'class'} (gradients of predicted classes w.r.t. inputs) or
#'    \code{'loss'} (gradients of loss w.r.t. inputs) (\code{'class'} by default).
#' @param top.n.genes Top n genes (positive and negative) taken per cell type. 
#'
#' @return Object containing gradients in the \code{interpret.gradients} slot of
#'   the \code{DeconvDLModel} object (\code{trained.model} slot).
#'
#' @export
#'
#' @seealso \code{\link{interGradientsDL}} \code{\link{trainDeconvModel}}
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(
#'     counts = matrix(
#'       rpois(30, lambda = 5), nrow = 15, ncol = 10,
#'       dimnames = list(paste0("Gene", seq(15)), paste0("RHC", seq(10)))
#'     )
#'   ),
#'   colData = data.frame(
#'     Cell_ID = paste0("RHC", seq(10)),
#'     Cell_Type = sample(x = paste0("CellType", seq(2)), size = 10,
#'                        replace = TRUE)
#'   ),
#'   rowData = data.frame(
#'     Gene_ID = paste0("Gene", seq(15))
#'   )
#' )
#' SDDLS <- createSpatialDDLSobject(
#'   sc.data = sce,
#'   sc.cell.ID.column = "Cell_ID",
#'   sc.gene.ID.column = "Gene_ID",
#' )
#' SDDLS <- genMixedCellProp(
#'   object = SDDLS,
#'   cell.ID.column = "Cell_ID",
#'   cell.type.column = "Cell_Type",
#'   num.sim.spots = 50,
#'   verbose = TRUE
#' )
#' SDDLS <- simMixedProfiles(SDDLS)
#' SDDLS <- trainDeconvModel(
#'   object = SDDLS,
#'   batch.size = 12,
#'   num.epochs = 5
#' )
#' ## calculating gradients
#' SDDLS <- interGradientsDL(SDDLS)
#' plotHeatmapGrads(SDDLS, top.n.genes = 2)
#' }
#'   
plotHeatmapGrads <- function(
    object, 
    method = "class", 
    top.n.genes = 15,
    scale.gradients = TRUE
  ) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) || 
      !requireNamespace("grid", quietly = TRUE)) {
    stop("ComplexHeatmap or grid R packages are required but not available")
  }
  if (!is(object, "SpatialDDLS")) {
    stop("The provided object is not of SpatialDDLS class")
  } else if (is.null(trained.model(object))) {
    stop("The provided object does not have a trained model for evaluation")
  }
  ## get metadata
  metadata.prop <- rbind(
    as.matrix(mixed.profiles(object, type.data = "train")@colData),
    as.matrix(mixed.profiles(object, type.data = "test")@colData)
  )
  metadata.prop <- metadata.prop[apply(X = metadata.prop == 100, MARGIN = 1, FUN = any), ]
  grads <- trained.model(object)@interpret.gradients[[method]]
  
  if (method == "class" | method == "loss") {
    top.genes <- topGradientsCellType(object, method = method, top.n.genes = top.n.genes)
    list.genes <- .sel.genes.sign(top.genes)
  } else {
    stop("method parameter has to be one of the following options: 'both', 'class' or 'loss'")
  }
  color.cell.types <- list(
    CellType = default.colors()[seq(ncol(metadata.prop))] %>% 
      setNames(colnames(metadata.prop))
  )
  named.vec <- sapply(
    X = seq(nrow(metadata.prop)),
    FUN = \(pos) {
      colnames(metadata.prop)[which(metadata.prop[pos, ] == 100)] %>% 
        setNames(rownames(metadata.prop)[pos])
    }
  ) 
  metadata.short <- data.frame(
    row.names = names(named.vec),
    CellType = named.vec
  )[rownames(grads), , drop = FALSE]
  ha <- ComplexHeatmap::HeatmapAnnotation(
    df = metadata.short,
    col = color.cell.types,
    annotation_name_gp = grid::gpar(fontsize = 10)
  )
  list.heatmaps <- list()
  for (i in c("Positive", "Negative")) {
    if (scale.gradients) {
      grads.f <- t(scale(grads[, list.genes[[i]]]))
    } else {
      grads.f <- t(grads[, list.genes[[i]]])
    }
    list.heatmaps[[i]] <- ComplexHeatmap::Heatmap(
      grads.f, 
      # col = greens,
      column_title = paste(i, "gradients (top", top.n.genes, "per cell type)"),
      top_annotation = ha, 
      border = "black", 
      row_names_gp = grid::gpar(fontsize = 8),
      name = "Score",
      column_title_gp = grid::gpar(fontface = "bold"),
      show_column_names = FALSE
    )
  }
  
  return(list.heatmaps)
}

.sel.genes.sign <- function(top.genes) {
  lapply(
    X = c("Positive", "Negative"), 
    FUN = \(sign.sel) {
      sapply(
        X = names(top.genes), 
        FUN = \(cell.type.sel) top.genes[[cell.type.sel]][[sign.sel]]
      ) %>% as.vector() %>% unique()
    }
  ) %>% setNames(c("Positive", "Negative"))
}

