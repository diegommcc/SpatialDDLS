################################################################################
############# Graph-based clustering using predicted proportions ###############
################################################################################

#' Cluster spatial data based on predicted cell proportions
#'
#' Cluster spatial transcriptomics data according to the cell proportions 
#' predicted in each spot. It allows to segregate data into niches with 
#' similar cell compositions. 
#'
#' @param object \code{\linkS4class{SpatialDDLS}} object with a trained 
#'    deconvolution model (\code{trained.model} slot) and pure mixed 
#'    transcriptional profiles (\code{mixed.profiles} slot).
#' @param index.st Name or index of the dataset/slide already deconvoluted to be
#'   clustered. If missing, all datasets already deconvoluted will be clustered.
#' @param method Clustering method. It can be \code{graph} (a nearest neighbor 
#'   graph is created and louvain algorithm is used to detect communities) or 
#'   \code{k.means} (k-means algorithm is run with the specified number of 
#'   centers (\code{k.centers} parameter)).
#' @param k.nn An integer specifying the number of nearest neighbors to during 
#'   graph construction (10 by default). Only if \code{method == "graph"}.
#' @param k.centers An integer specifying the number of centers for k-means 
#' algorithm. Only if \code{method == "k.means"}.
#' @param verbose Show informative messages during the execution (\code{TRUE} by
#'   default).
#'
#' @return A \code{\linkS4class{SpatialDDLS}} object containing calculated 
#'   clusters as a column in the slot \code{colData} of the 
#'   \code{\linkS4class{SpatialExperiment}} objects. 
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
#' ) %>% genMixedCellProp(
#'   object = SDDLS,
#'   cell.ID.column = "Cell_ID",
#'   cell.type.column = "Cell_Type",
#'   num.sim.spots = 50,
#'   verbose = TRUE
#' ) %>%simMixedProfiles(SDDLS) %>% 
#' trainDeconvModel(
#'   object = SDDLS,
#'   batch.size = 12,
#'   num.epochs = 5
#' ) %>% deconvSpatialDDLS(
#'   object = SDDLS,
#'   index.st = 1
#' ) %>% spatialPropClustering(index.st = 1)
#' }
#'   
spatialPropClustering <- function(
    object,
    index.st,
    method = "graph",
    k.nn = 10,
    k.centers = 5,
    verbose = TRUE
) {
  if (!requireNamespace("bluster", quietly = TRUE)) {
    stop("bluster R package is required but not available")
  }
  if (!is(object, "SpatialDDLS")) {
    stop("The provided object is not of SpatialDDLS class")
  } else if (is.null(spatial.experiments(object))) {
    stop("`spatial.experiments` slot is empty")
  } else if (is.null(deconv.spots(object))) {
    stop("`desconv.slots` slot is empty")
  }
  ## checking index.st
  if (missing(index.st)) {
    if (verbose) {
      message(
        "   No 'index.st' provided. Deconvoluting all SpatialExperiment objects ", 
        "contained in the `spatial.experiments` slot\n"
      ) 
    }
    index.st <- seq_along(spatial.experiments(object))
  } else {
    if (is.character(index.st) & !is.null(names(spatial.experiments(object)))) {
      ## check if all index.st are present in the slot
      stopifnot(
        "`index.st` contains elements not present in spatial.experiments slot " = index.st %in% 
          names(spatial.experiments(object))
      )
    }
  }
  ## I should check this: when objects are named
  # if (is.character(index.st)) {
  #   namesList <- index.st
  # } else {
  #   namesList <- names(spatial.experiments(object))[index.st]
  # }
  ## chceking predicted cell proportions
  if (is.null(deconv.spots(object, index.st = index.st))) {
    stop("No predicted cell proportions were provided for the selected ST dataset")
  } else if (is.null(spatial.experiments(object, index.st = index.st))) { ## not sure if this is needed
    stop("No spatial data")
  }
  ## clustering method
  if (method == "graph") {
    method.clustering <- bluster::NNGraphParam(k = k.nn, cluster.fun = "louvain")
    if (verbose) message("=== Selected graph-based clustering\n")
  } else if (method == "k.means") {
    method.clustering <- bluster::KmeansParam(centers = k.centers)
    if (verbose) message("=== Selected k-means clustering\n")
  } else {
    stop("method for clustering not available. Possible options are 'graph' and 'k.means'")
  }
  ## clustering: it has to be vectorized to take into account several slides
  clustering.res <- lapply(
    index.st, 
    FUN = \(index.pointer) {
      if (verbose) message(paste("=== Running clustering for slide", index.pointer))
      res <- bluster::clusterRows(
        x = deconv.spots(object, index.st = index.pointer)[["Regularized"]], 
        BLUSPARAM = method.clustering
      )    
      return(factor(paste0("C", res), levels = paste0("C", levels(res))))
    }
  )
  ## introducing results into the Spatial objects
  spatial.experiments(object) <- lapply(
    X = index.st,
    FUN = \(st.pointer) {
      st.obj <- spatial.experiments(object, index.st = st.pointer)
      colData(st.obj)[[paste(
        "Clustering", method, "k",
        ifelse(method == "graph", k.nn, k.centers), sep = "."
      )]]  <- clustering.res[[st.pointer]]
      return(st.obj)
    }
  ) 
  
  return(object)
}
