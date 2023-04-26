#' @importFrom ggplot2 scale_color_gradientn scale_fill_gradientn
NULL
color.prop.scale.blues <- c(
  "#ECF4FB", "#E1EDF8", "#D7E6F4", "#CDE0F1", "#C1D9ED", "#B0D2E7", "#A0CAE1",
  "#8BBFDC", "#75B3D8", "#62A8D2", "#519CCB", "#4090C5", "#3282BD", "#2474B6", 
  "#1966AD", "#0E59A2", "#084B94", "#083D7F", "#08306B"
)
color.prop.scale.spectral <- grDevices::colorRampPalette(
  colors = rev(x = RColorBrewer::brewer.pal(n = 11, name = "Spectral"))
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
#'   can be either a position or a name if a named list was provided.
#' @param colors Color scale to be used. It can be \code{"blues"} or
#'   \code{"spectral"} (the former by default).
#' @param set If results were simplified (see \code{?\link{deconvSpatialDDLS}}
#'   for details), what results to plot (\code{raw} by default).
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
  if(is.list(st.pred)) {
    if (!set %in% c("raw", "simplify.set", "simpli.majority")) {
      stop(
        "`set`must be one of the following options: 'raw', 'simplify.set', 'simpli.majority'"
      )
    }
    st.pred <- st.pred[[set]]
  }
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
#'   can be either a position or a name if a named list was provided.
#' @param cell.type Cell type predicted proportions to color spots by.
#' @param colors Color scale to be used. It can be \code{"blues"} or
#'   \code{"spectral"} (the former by default).
#' @param set If results were simplified (see \code{?\link{deconvSpatialDDLS}}
#'   for details), what results to plot (\code{raw} by default).
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
  if(is.list(st.pred)) {
    if (!set %in% c("raw", "simplify.set", "simpli.majority")) {
      stop(
        "`set`must be one of the following options: 'raw', 'simplify.set', 'simpli.majority'"
      )
    }
    st.pred <- st.pred[[set]]
  }
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
    
  if (is.null(title)) title.plot <- paste0("Predicted proportions (", cell.type, ")")
  
  plot <- ggplot(
    dfPlot, aes(
      x = .data[["Spatial 1"]], y = .data[["Spatial 2"]], 
      color = .data[[cell.type]]
    )
  ) + geom_point(size = size.point) + scale_colors + 
    ggtitle(title.plot) + SpatialDDLSTheme() 

  return(plot)
}

