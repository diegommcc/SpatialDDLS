

color.prop.scale <- c(
  "#ECF4FB", "#E1EDF8", "#D7E6F4", "#CDE0F1", "#C1D9ED", "#B0D2E7", "#A0CAE1",
  "#8BBFDC", "#75B3D8", "#62A8D2", "#519CCB", "#4090C5", "#3282BD", "#2474B6", 
  "#1966AD", "#0E59A2", "#084B94", "#083D7F", "#08306B"
)

################################################################################
####################### Plot spatial proportions (all) #########################
################################################################################

#' Plot predicted proportions for all cell types using spatial coordinates of
#' spots
#'
#' Color spots on the spatial coordinates plot according to their predicted cell
#' type proportions. All cell types are represented altogether using the same
#' scale of colors (from 0 to 1).
#'
#' @param object \code{\linkS4class{SpatialDDLS}} object.
#' @param index.st
#' @param colors Vector of colors to be used. Only vectors with a number of
#'   colors equal to or greater than the levels of \code{color.by} will be
#'   accepted. By default, a custom color list is used.
#' @param set
#' @param color.scale Variable used to display data in different panels. If
#'   \code{NULL}, the plot is not split into different panels. Options are
#'   \code{nCellTypes} (by number of different cell types) and \code{CellType}
#'   (by cell type).
#' @param size.point Size of points (0.1 by default).
#' @param alpha.point Alpha of points (0.1 by default).
#' @param title Boolean indicating whether single-cell profiles are filtered out
#'   and only errors associated with pseudo-bulk samples are displayed
#'   (\code{TRUE} by default).
#' @param nrow Number of rows in the split plot.
#' @param ncol Nomber of columns in the split plot
#' @param theme \pkg{ggplot2} theme.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @seealso \code{\link{plotSpatialProp}} \code{\link{deconvSpatialDDLS}}
#'   \code{\link{trainDeconvModel}} 
#'
#' @examples
#' \dontrun{
#'
#' }
#' 
plotSpatialPropAll <- function(
    object,
    index.st,
    colors = NULL,
    set = "raw",
    size.point = 0.1,
    alpha.point = 1,
    title = NULL,
    nrow = NULL,
    ncol = NULL,
    theme = NULL,
    ...
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
  ## TODO: change the default behaviour: create a list with diff elements
  ## also consider the possibitly of using raw, simplified props, etc
  st.pred <- deconv.spots(object = object, index.st = index.st)
  if(is.list(st.pred)) {
    if (!set %in% c("raw", "simplify.set", "simpli.majority")) {
      stop(
        "`set`must be one of the following options: 'raw', 'simplify.set', 'simpli.majority'"
      )
    }
    st.pred <- st.pred[[set]]
  }
  # if (!all(rownames(st.coor) == rownames(st.pred))) {
  #   stop("Matrices in `deconv.spots` and `spatalCoords` are not complementary")
  # }
  dfPlot <- reshape2::melt(
    as.data.frame(cbind(st.coor, st.pred)), 
    id.vars = c("Spatial 1", "Spatial 2"), 
    variable.name = "CellType", value.name = "Proportion"
  )
  # TODO: change this to parse a gradient palette blablabla
  if (is.null(colors))  {
    scale_colors <- scale_color_gradientn(colors = color.prop.scale, limits = c(0, 1))  
    # colors <- scale_color_gradient(low = "white", high = "blue")
  } else {
    scale_colors <- scale_color_gradientn(colors = colors, limits = c(0, 1))  
  }
    
  if (is.null(title)) title <- "Predicted proportions"
  plot <- ggplot(
    dfPlot, aes(x = .data[["Spatial 1"]], y = .data[["Spatial 2"]], color = Proportion)
  ) + geom_point(size = size.point, alpha = alpha.point) + 
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
#' Color spots on the spatial coordinates plot according to predicted
#' proportions of a particular cell type. Color scale is adapted depending on
#' the range of predicted proportions.
#'
#' @param object \code{\linkS4class{SpatialDDLS}} object.
#' @param index.st Nmae or index of the desired
#'   \code{\linkS4class{SpatialExperiment}} object to visualize.
#' @param cell.type Cell type predicted proportions to color spots by.
#' @param set
#' @param color.scale
#' @param size.point Size of points (0.1 by default).
#' @param alpha.point Alpha of points (0.1 by default).
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
#' @examples
#' \dontrun{
#'
#' }
#' 
plotSpatialProp <- function(
    object,
    index.st,
    cell.type,
    colors = NULL, 
    set = "raw",
    size.point = 1,
    alpha.point = 1,
    title = NULL,
    theme = NULL,
    ...
) {
  if (!is(object, "SpatialDDLS")) {
    stop("The provided object is not of class SpatialDDLS")
  } else if (is.null(spatial.experiments(object)) || is.null(deconv.spots(object))) {
    stop("djhdhdhd")
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
  
  # if (!all(rownames(st.coor) == rownames(st.pred))) {
  #   stop("Matrices in `deconv.spots` and `spatalCoords` are not complementary")
  # }
  if (is.null(colors)) {
    color.scale <- scale_color_gradientn(colors = color.prop.scale)
    # color.scale <- scale_color_gradient(low = "white", high = "blue")
  } else {
    color.scale <- scale_color_gradientn(colors = colors)
  }
    
  if (is.null(title)) title.plot <- paste0("Predicted proportions (", cell.type, ")")
  
  plot <- ggplot(
    dfPlot, aes(
      x = .data[["Spatial 1"]], y = .data[["Spatial 2"]], 
      color = .data[[cell.type]]
    )
  ) + geom_point(size = size.point, alpha = alpha.point) + color.scale + 
    ggtitle(title.plot) + SpatialDDLSTheme() 

  return(plot)
}

################################################################################
##################### Plot spatial proportions (blended) #######################
################################################################################

# not exported, still unstable
.plotSpatialPropBlended <- function(
    object,
    index.st,
    colors,
    cell.type,
    set = "raw",
    size.point = 2,
    alpha.point = 0.5,
    title = NULL,
    theme = NULL,
    ...
) {
  if (!is(object, "SpatialDDLS")) {
    stop("The provided object is not of class SpatialDDLS")
  } else if (is.null(spatial.experiments(object)) || is.null(deconv.spots(object))) {
    stop("djhdhdhd")
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
  # if (!all(rownames(st.coor) == rownames(st.pred))) {
  #   stop("Matrices in `deconv.spots` and `spatalCoords` are not complementary")
  # }
  if (length(cell.type) != 2) stop("Only 2 cell types are accepted")
  dfPlot <- as.data.frame(cbind(st.coor, st.pred[, cell.type]))
  if (is.null(title)) title.plot <- "Predicted proportions"
  
  plot <- ggplot(
    dfPlot, aes(x = .data[["Spatial 1"]], y = .data[["Spatial 2"]])
  ) + geom_point(
    aes(fill = .data[[cell.type[1]]], color = .data[[cell.type[1]]]), 
    size = size.point, alpha = alpha.point
  ) + scale_color_gradientn(cell.type[1], colors = c("white", "blue")) + 
    scale_fill_gradientn(cell.type[1], colors = c("white", "blue")) + 
    ggnewscale::new_scale_color() +
    geom_point(
      aes(fill = .data[[cell.type[2]]]),  color = "black",
      size = size.point, alpha = alpha.point, shape = 21
    ) + 
    scale_fill_gradientn(cell.type[2], colors = c("white", "red")) + 
    SpatialDDLSTheme() 
  
  return(plot)
}
