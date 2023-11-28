#' @importFrom ggplot2 scale_color_gradientn scale_color_gradient2
#' @importFrom stats median
NULL
color.prop.scale.blues <- c(
  "#ECF4FB", "#E1EDF8", "#D7E6F4", "#CDE0F1", "#C1D9ED", "#B0D2E7", "#A0CAE1",
  "#8BBFDC", "#75B3D8", "#62A8D2", "#519CCB", "#4090C5", "#3282BD", "#2474B6", 
  "#1966AD", "#0E59A2", "#084B94", "#083D7F", "#08306B"
)
color.prop.scale.spectral <- grDevices::colorRampPalette(
  colors = rev(
    c(
      "#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", 
      "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2"
    )
  )
)(100)

################################################################################
####################### Plot spatial proportions (all) #########################
################################################################################

#' Plot predicted proportions for all cell types using spatial coordinates of
#' spots
#'
#' Color spots on the spatial coordinates plot according to their predicted cell
#' type proportions. All cell types are represented together using the same
#' color scale from 0 to 1.
#'
#' @param object A \code{\linkS4class{SpatialDDLS}} object.
#' @param index.st Index of the spatial transcriptomics data to be plotted. It
#'   can be either a position or a name if a named list of 
#'   \code{\linkS4class{SpatialExperiment}} objects was provided.
#' @param colors Color scale to be used. It can be \code{"blues"} or
#'   \code{"spectral"} (the former by default).
#' @param set If results were simplified (see \code{?\link{deconvSpatialDDLS}}
#'   for details), which results to plot (\code{raw} by default).
#' @param prediction It can be \code{"Regularized"}, \code{"Intrinsic"} or 
#'   \code{"Extrinsic"} (\code{"Regularized"} by default).
#' @param size.point Size of points (0.1 by default).
#' @param title Title of plot.
#' @param nrow Number of rows in the split plot.
#' @param ncol Number of columns in the split plot.
#' @param theme \pkg{ggplot2} theme.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @seealso \code{\link{plotSpatialProp}} \code{\link{deconvSpatialDDLS}}
#'   \code{\link{trainDeconvModel}}
#' 
plotSpatialPropAll <- function(
    object,
    index.st,
    colors = "blues",
    set = "raw",
    prediction = "Regularized",
    size.point = 0.1,
    title = NULL,
    nrow = NULL,
    ncol = NULL,
    theme = NULL
) {
  if (!is(object, "SpatialDDLS")) {
    stop("The provided object is not of class SpatialDDLS")
  } else if (
    is.null(spatial.experiments(object)) || is.null(deconv.spots(object))
  ) {
    stop(
      "Either predictions (`deconv.spots` slot) or ST data ", 
      "(`spatial.experiments` slot) not present in SpatialDDLS object"
    )
  } 
  
  ## getting data
  st.coor <- SpatialExperiment::spatialCoords(
    spatial.experiments(object = object, index.st = index.st)
  )[, 1:2]
  colnames(st.coor) <- paste("Spatial", 1:2)
  ## TODO: change the default behavior: create a list with diff elements
  ## also consider the possibility of using raw, simplified props, etc
  st.pred <- deconv.spots(object = object, index.st = index.st)
  if (is.list(st.pred) & any(names(st.pred) %in% c("raw", "simplify.set", "simpli.majority"))) {
    if (!set %in% c("raw", "simplify.set", "simpli.majority")) {
      stop(
        "`set`must be one of the following options: 'raw', 'simplify.set', 'simpli.majority'"
      )
    }
    st.pred <- st.pred[[set]]
  }
  if (!any(prediction %in% c("Regularized", "Intrinsic", "Extrinsic"))) {
    stop("prediction can only be one of the following options: 'Regularized', 'Intrinsic', 'Extrinsic'")
  }
  st.pred <- st.pred[[prediction]]
  
  dfPlot <- reshape2::melt(
    as.data.frame(cbind(st.coor, st.pred)), 
    id.vars = c("Spatial 1", "Spatial 2"), 
    variable.name = "CellType", value.name = "Proportion"
  )
  if (colors == "blues")  {
    scale_colors <- scale_color_gradientn(
      colors = color.prop.scale.blues, limits = c(0, 1)
    )  
    # colors <- scale_color_gradient(low = "white", high = "blue")
  } else if (colors == "spectral") {
    scale_colors <- scale_color_gradientn(
      colors = color.prop.scale.spectral, limits = c(0, 1)
    )  
  } else {
    scale_colors <- scale_color_gradientn(
      colors = color.prop.scale.blues, limits = c(0, 1)
    )
  }
    
  if (is.null(title)) title <- "Predicted proportions"
  plot <- ggplot(
    dfPlot, 
    aes(
      x = .data[["Spatial 1"]], y = .data[["Spatial 2"]], 
      color = .data[["Proportion"]]
    )
  ) + geom_point(size = size.point) + 
    ggtitle(title) + 
    SpatialDDLSTheme() + facet_wrap(~ CellType, nrow = nrow, ncol = ncol) +
    scale_colors
  
  return(plot)
}


################################################################################
##################### Plot spatial proportions (single) ########################
################################################################################

#' Plot predicted proportions for a specific cell type using spatial coordinates
#' of spots
#'
#' Color spots on the spatial coordinates according to the predicted proportions
#' of a particular cell type. Color scale is adapted depending on the range of
#' predicted proportions.
#'
#' @param object A \code{\linkS4class{SpatialDDLS}} object.
#' @param index.st Index of the spatial transcriptomics data to be plotted. It
#'   can be either a position or a name if a named list of 
#'   \code{\linkS4class{SpatialExperiment}} objects was provided.
#' @param cell.type Cell type predicted proportions to color spots by.
#' @param colors Color scale to be used. It can be \code{"blues"} or
#'   \code{"spectral"} (the former by default).
#' @param set If results were simplified (see \code{?\link{deconvSpatialDDLS}}
#'   for details), what results to plot (\code{raw} by default).
#' @param prediction It can be \code{"Regularized"}, \code{"Intrinsic"} or 
#'   \code{"Extrinsic"} (\code{"Regularized"} by default).
#' @param limits A vector of two elements indicating wanted limits for color
#'   scale. If \code{NULL} (by default), color scale is adjusted to max and min
#'   predicted proportions.
#' @param size.point Size of points (0.1 by default).
#' @param title Title of plot.
#' @param theme \pkg{ggplot2} theme.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @seealso \code{\link{plotSpatialPropAll}} \code{\link{deconvSpatialDDLS}}
#'   \code{\link{trainDeconvModel}}
#'   
plotSpatialProp <- function(
    object,
    index.st,
    cell.type,
    colors = "blues", 
    set = "raw",
    prediction = "Regularized",
    limits = NULL,
    size.point = 1,
    title = NULL,
    theme = NULL
) {
  if (!is(object, "SpatialDDLS")) {
    stop("The provided object is not of class SpatialDDLS")
  } else if (is.null(spatial.experiments(object)) || is.null(deconv.spots(object))) {
    stop(
      "Either predictions (`deconv.spots` slot) or ST data ", 
      "(`spatial.experiments` slot) not present in SpatialDDLS object"
    )
  } 
  ## getting data
  st.coor <- SpatialExperiment::spatialCoords(
    spatial.experiments(object = object, index.st = index.st)
  )[, 1:2]
  colnames(st.coor) <- paste("Spatial", 1:2)
  st.pred <- deconv.spots(object = object, index.st = index.st)
  if(is.list(st.pred) & any(names(st.pred) %in% c("raw", "simplify.set", "simpli.majority"))) {
    if (!set %in% c("raw", "simplify.set", "simpli.majority")) {
      stop(
        "`set`must be one of the following options: 'raw', 'simplify.set', 'simpli.majority'"
      )
    }
    st.pred <- st.pred[[set]]
  }
  if (!any(prediction %in% c("Regularized", "Intrinsic", "Extrinsic"))) {
    stop("prediction can only be one of the following options: 'Regularized', 'Intrinsic', 'Extrinsic'")
  }
  st.pred <- st.pred[[prediction]]
  
  
  if (!cell.type %in% colnames(st.pred)) stop("`cell.type` must be a valid cell type")
  
  st.pred <- st.pred[, cell.type, drop = FALSE]
  dfPlot <- as.data.frame(cbind(st.coor, st.pred))
  if (colors == "blues")  {
    scale_colors <- scale_color_gradientn(
      colors = color.prop.scale.blues, limits = limits
    )  
    # colors <- scale_color_gradient(low = "white", high = "blue")
  } else if (colors == "spectral") {
    scale_colors <- scale_color_gradientn(
      colors = color.prop.scale.spectral, limits = limits
    )  
  } else {
    scale_colors <- scale_color_gradientn(
      colors = color.prop.scale.blues, limits = limits
    )
  }
    
  if (is.null(title)) title <- paste0("Predicted proportions (", cell.type, ")")
  
  plot <- ggplot(
    dfPlot, aes(
      x = .data[["Spatial 1"]], y = .data[["Spatial 2"]], 
      color = .data[[cell.type]]
    )
  ) + geom_point(size = size.point) + scale_colors + 
    ggtitle(title) + SpatialDDLSTheme() 

  return(plot)
}

################################################################################
########################### Plot gene expression ###############################
################################################################################

#' Plot normalized gene expression data (logCPM) in spatial coordinates
#'
#' Color spots on the spatial coordinates according to the logCPM values of a 
#' particular gene.
#'
#' @param object A \code{\linkS4class{SpatialDDLS}} object.
#' @param index.st Index of the spatial transcriptomics data to be plotted. It
#'   can be either a position or a name if a named list of 
#'   \code{\linkS4class{SpatialExperiment}} objects was provided.
#' @param gene Gene to color spots by.
#' @param colors Color scale to be used. It can be \code{"blues"} or
#'   \code{"spectral"} (the latter by default).
#' @param size.point Size of points (0.1 by default).
#' @param title Title of plot.
#' @param theme \pkg{ggplot2} theme.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @seealso \code{\link{interGradientsDL}} \code{\link{topGradientsCellType}}
#'   
plotSpatialGeneExpr <- function(
    object,
    index.st,
    gene,
    colors = "spectral", 
    size.point = 1,
    title = NULL,
    theme = NULL
) {
  if (!is(object, "SpatialDDLS")) {
    stop("The provided object is not of class SpatialDDLS")
  } 
  if (missing(index.st)) {
    message(
      "   No 'index.st' provided. Using first ST dataset"
    ) 
    index.st <- 1
  } else {
    if (is.character(index.st) & !is.null(names(spatial.experiments(object)))) {
      ## check if all index.st are present in the slot
      stopifnot(
        "`index.st` contains elements not present in spatial.experiments slot " = index.st %in% 
          names(spatial.experiments(object))
      )
    }
  }
  ## getting data
  st.coor <- SpatialExperiment::spatialCoords(
    spatial.experiments(object = object, index.st = index.st)
  )[, 1:2]
  colnames(st.coor) <- paste("Spatial", 1:2)
  ## checking if gene is in data
  if (!any(gene %in% rownames(assays(
    spatial.experiments(object = object, index.st = index.st)
  )[[1]]))) {
    stop("Provided gene is not in spatial data")
  }
  gene.expr <- assays(
    spatial.experiments(object = object, index.st = index.st)
  )[[1]][gene, , drop = FALSE] %>% as.matrix()
  
  gene.expr.norm <- log2(.cpmCalculate(x = t(gene.expr) + 1))
  
  gene.expr.norm <- gene.expr.norm[rownames(st.coor), , drop = FALSE]
  dfPlot <- as.data.frame(cbind(st.coor, gene.expr.norm))
  if (colors == "blues")  {
    scale_colors <- scale_color_gradientn(
      colors = color.prop.scale.blues
    )  
    # colors <- scale_color_gradient(low = "white", high = "blue")
  } else if (colors == "spectral") {
    scale_colors <- scale_color_gradientn(
      colors = color.prop.scale.spectral
    )  
  } else {
    scale_colors <- scale_color_gradientn(
      colors = color.prop.scale.blues
    )
  }
  
  if (is.null(title)) title <- paste0("Normalized gene expression (", gene, ")")
  
  plot <- ggplot(
    dfPlot, aes(
      x = .data[["Spatial 1"]], y = .data[["Spatial 2"]], 
      color = .data[[gene]]
    )
  ) + geom_point(size = size.point) + scale_colors + 
    ggtitle(title) + SpatialDDLSTheme() 
  
  return(plot)
}

################################################################################
###################### Plot spatial clusters (single) ##########################
################################################################################

#' Plot results of clustering based on predicted cell proportions
#'
#' Color spots on the spatial coordinates according to the results of 
#' clustering based on predicted proportions. 
#'
#' @param object A \code{\linkS4class{SpatialDDLS}} object.
#' @param index.st Index of the spatial transcriptomics data to be plotted. It
#'   can be either a position or a name if a named list of 
#'   \code{\linkS4class{SpatialExperiment}} objects was provided.
#' @param method Clustering method results to plot. It can be \code{"graph"} or
#'   \code{"k.means"}. If missing, the first configuration found in the object 
#'   will be plotted. 
#' @param k.nn Number of nearest neighbors used if \code{ method == "graph"}. 
#' @param k.centers Number of k centers used if \code{ method == "k.means"}. 
#' @param colors Vector of colors to be used.
#' @param size.point Size of points (0.1 by default).
#' @param title Title of plot.
#' @param theme \pkg{ggplot2} theme.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @seealso \code{\link{spatialPropClustering}} \code{\link{deconvSpatialDDLS}}
#'   
plotSpatialClustering <- function(
    object,
    index.st,
    method,
    k.nn,
    k.centers,
    colors, 
    size.point = 1,
    title = NULL,
    theme = NULL
) {
  if (!is(object, "SpatialDDLS")) {
    stop("The provided object is not of class SpatialDDLS")
  } else if (is.null(spatial.experiments(object))) {
    stop("`spatial.experiments` slot is empty")
  }
  ## checking index
  if (missing(index.st)) {
    message(
      "   No 'index.st' provided. Using first ST dataset"
    ) 
    index.st <- 1
  } else {
    if (is.character(index.st) & !is.null(names(spatial.experiments(object)))) {
      ## check if all index.st are present in the slot
      stopifnot(
        "`index.st` contains elements not present in spatial.experiments slot " = index.st %in% 
          names(spatial.experiments(object))
      )
    }
  }
  ## getting data
  st.coor <- SpatialExperiment::spatialCoords(
    spatial.experiments(object = object, index.st = index.st)
  )[, 1:2]
  colnames(st.coor) <- paste("Spatial", 1:2)
  st.clusternig <- SummarizedExperiment::colData(
    spatial.experiments(object = object, index.st = index.st)
  )
  if (missing(method)) {
    cls <- grep(
      pattern = "^Clustering\\.", x = colnames(st.clusternig), value = TRUE
    )[1]
    message(paste0("=== Plotting first clustering configuration ", cls, "\n"))
  } else if (method == "graph") {
    if (missing(k.nn)) {
      cls <- grep(
        pattern = "^Clustering\\.graph\\.", 
        x = colnames(st.clusternig), value = TRUE
      )[1]
      message(
        paste0("=== Plotting first graph clustering configuration ", cls, "\n")
      )
    } else {
      cls <- paste0("Clustering.graph.k.", k.nn)
    }
  } else if (method == "k.means") {
    if (missing(k.centers)) {
      cls <- grep(
        pattern = "^Clustering\\.k.means\\.", x = colnames(st.clusternig), 
        value = TRUE
      )[1]
      message(
        paste0("=== Plotting first k-means clustering configuration ", cls, "\n")
      )
    } else {
      cls <- paste0("Clustering.k.means.k.", k.centers)
    }
  }
  ## check if selected config exists
  if (is.null(st.clusternig[[cls]])) {
    stop("Selected clustering configuration does not exist")
  }
  st.pred <- st.clusternig[, cls, drop = FALSE]
  dfPlot <- as.data.frame(cbind(st.coor, st.pred))
  ## colors
  if (missing(colors)) colors <- default.colors()
  if (length(colors) < length(unique(st.pred[[cls]]))) 
    stop("Number of provided colors is not large enough")
  
  if (is.null(title)) title <- paste0("Clustering results (", cls, ")")
  
  plot <- ggplot(
    dfPlot, aes(
      x = .data[["Spatial.1"]], y = .data[["Spatial.2"]], 
      color = .data[[cls]]
    )
  ) + geom_point(size = size.point) + 
    scale_color_manual(values = colors, name = cls) + 
    ggtitle(title) + SpatialDDLSTheme() 
  
  return(plot)
}

################################################################################
################################ Plot distances ################################
################################################################################

#' Plot distances between intrinsic and extrinsic profiles
#'
#' Color spots on the spatial coordinates according to distances between 
#' intrinsic and extrinsic transcriptional profiles. 
#'
#' @param object A \code{\linkS4class{SpatialDDLS}} object.
#' @param index.st Index of the spatial transcriptomics data to be plotted. It
#'   can be either a position or a name if a named list was provided.
#' @param mid.scale The midpoint of the diverging scale. it may be \code{'mean'}
#'   or \code{'median'} (the former by default).
#' @param size.point Size of points (0.1 by default).
#' @param title Title of plot.
#' @param theme \pkg{ggplot2} theme.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @seealso \code{\link{deconvSpatialDDLS}} \code{\link{trainDeconvModel}}
#'
plotDistances <- function(
    object,
    index.st,
    mid.scale = "mean",
    size.point = 1,
    title = NULL,
    theme = NULL
) {
  if (!is(object, "SpatialDDLS")) {
    stop("The provided object is not of class SpatialDDLS")
  } else if (is.null(spatial.experiments(object))) {
    stop("`spatial.experiments` slot is empty")
  }
  ## checking index
  if (missing(index.st)) {
    message(
      "   No 'index.st' provided. Using first ST dataset"
    )
    index.st <- 1
  } else {
    if (is.character(index.st) & !is.null(names(spatial.experiments(object)))) {
      stopifnot(
        "`index.st` contains elements not present in spatial.experiments slot " = index.st %in%
          names(spatial.experiments(object))
      )
    }
  }
  ## getting data
  st.coor <- SpatialExperiment::spatialCoords(
    spatial.experiments(object = object, index.st = index.st)
  )[, 1:2]
  colnames(st.coor) <- paste("Spatial", 1:2)
  st.distances <- deconv.spots(object, index.st = index.st)[["Distances"]]
  dfPlot <- as.data.frame(cbind(st.coor, Distances = st.distances))
  ## colors
  if (is.null(title)) 
    title <- paste0("Distances between intrinsic and extrinsic profiles")
  
  if (mid.scale == "mean") {
    midpoint <- mean(dfPlot[["Distances"]])
  } else if (mid.scale == "median") {
    midpoint <- median(dfPlot[["Distances"]])
  } else {
    # warning("Using mean as mean distances as mid point")
    midpoint <- mean(dfPlot[["Distances"]])
  }
  plot <- ggplot(
    dfPlot, aes(
      x = .data[["Spatial 1"]], y = .data[["Spatial 2"]],
      color = .data[["Distances"]]
    )
  ) + geom_point(size = size.point) +
    scale_color_gradient2(
      midpoint =  midpoint, 
      low = "blue", mid = "white",
      high = "red", space = "Lab" 
    ) + ggtitle(title) + SpatialDDLSTheme()

  return(plot)
}
