

################################################################################
####################### Plot spatial proportions (all) #########################
################################################################################

#' Generate correlation plots between predicted and expected cell type
#' proportions from test data
#'
#' Generate correlation plot between predicted and expected cell type
#' proportions from test data. Correlation plots can be displayed all mixed or
#' split by cell type (\code{CellType}) or number of cell types present in the
#' samples (\code{nCellTypes}). See the \code{facet.by} argument and examples
#' for more information. Moreover, a user-selected correlation value is
#' displayed as an annotation on the plots. See the \code{corr} argument for
#' details.
#'
#' @param object \code{\linkS4class{SpatialDDLS}} object with
#'   \code{trained.model} slot containing metrics in the
#'   \code{test.deconv.metrics} slot of a
#'   \code{\linkS4class{SpatialDDLSDNN}} object.
#' @param colors Vector of colors to be used. Only vectors with a number of
#'   colors equal to or greater than the levels of \code{color.by} will be
#'   accepted. By default, a custom color list is used.
#' @param facet.by Variable used to display data in different panels. If
#'   \code{NULL}, the plot is not split into different panels. Options are
#'   \code{nCellTypes} (by number of different cell types) and \code{CellType}
#'   (by cell type).
#' @param color.by Variable used to color data. Options are \code{nCellTypes}
#'   and \code{CellType}.
#' @param corr Correlation value displayed as an annotation on the plot.
#'   Available metrics are Pearson's correlation coefficient (\code{'pearson'})
#'   and concordance correlation coefficient (\code{'ccc'}). The argument can be
#'   \code{'pearson'}, \code{'ccc'} or \code{'both'} (by default).
#' @param filter.sc Boolean indicating whether single-cell profiles are filtered
#'   out and only errors associated with pseudo-bulk samples are displayed
#'   (\code{TRUE} by default).
#' @param pos.x.label X-axis position of correlation annotations (0.95 by
#'   default).
#' @param pos.y.label Y-axis position of correlation annotations (0.1 by
#'   default).
#' @param sep.labels Space separating annotations if \code{corr} is equal to
#'   \code{'both'} (0.15 by default).
#' @param size.point Size of points (0.1 by default).
#' @param alpha.point Alpha of points (0.1 by default).
#' @param nrow Number of rows if \code{facet.by} is different from \code{NULL}.
#' @param ncol Number of columns if \code{facet.by} is other than \code{NULL}.
#' @param title Title of the plot.
#' @param theme \pkg{ggplot2} theme.
#' @param ... Additional arguments for the \link[ggplot2]{facet_wrap} function
#'   from \pkg{ggplot2} if \code{facet.by} is not \code{NULL}.
#'
#' @return A ggplot object with the correlation plots between expected and
#'   actual proportions.
#'
#' @export
#'
#' @seealso \code{\link{calculateEvalMetrics}} \code{\link{distErrorPlot}}
#'   \code{\link{blandAltmanLehPlot}} \code{\link{barErrorPlot}}
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(
#'     counts = matrix(
#'       rpois(30, lambda = 5), nrow = 15, ncol = 20,
#'       dimnames = list(paste0("Gene", seq(15)), paste0("RHC", seq(20)))
#'     )
#'   ),
#'   colData = data.frame(
#'     Cell_ID = paste0("RHC", seq(20)),
#'     Cell_Type = sample(x = paste0("CellType", seq(6)), size = 20,
#'                        replace = TRUE)
#'   ),
#'   rowData = data.frame(
#'     Gene_ID = paste0("Gene", seq(15))
#'   )
#' )
#' DDLS <- loadSCProfiles(
#'   single.cell.data = sce,
#'   cell.ID.column = "Cell_ID",
#'   gene.ID.column = "Gene_ID"
#' )
#' probMatrixValid <- data.frame(
#'   Cell_Type = paste0("CellType", seq(6)),
#'   from = c(1, 1, 1, 15, 15, 30),
#'   to = c(15, 15, 30, 50, 50, 70)
#' )
#' DDLS <- generateBulkCellMatrix(
#'   object = DDLS,
#'   cell.ID.column = "Cell_ID",
#'   cell.type.column = "Cell_Type",
#'   prob.design = probMatrixValid,
#'   num.bulk.samples = 50,
#'   verbose = TRUE
#' )
#' # training of DDLS model
#' tensorflow::tf$compat$v1$disable_eager_execution()
#' DDLS <- trainSpatialDDLSModel(
#'   object = DDLS,
#'   on.the.fly = TRUE,
#'   batch.size = 15,
#'   num.epochs = 5
#' )
#' # evaluation using test data
#' DDLS <- calculateEvalMetrics(
#'   object = DDLS
#' )
#' # correlations by cell type
#' corrExpPredPlot(
#'   object = DDLS,
#'   facet.by = "CellType",
#'   color.by = "CellType",
#'   corr = "both"
#' )
#' # correlations of all samples mixed
#' corrExpPredPlot(
#'   object = DDLS,
#'   facet.by = NULL,
#'   color.by = "CellType",
#'   corr = "ccc",
#'   pos.x.label = 0.2,
#'   alpha.point = 0.3
#' )
#' }
#' 
plotSpatialPropAll <- function(
    object,
    index.st,
    colors,
    set = "raw",
    color.scale = "shared",
    size.point = 0.1,
    alpha.point = 1,
    title = NULL,
    theme = NULL,
    nrow = NULL,
    ncol = NULL,
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
  if (!all(rownames(st.coor) == rownames(st.pred))) {
    stop("Matrices in `deconv.spots` and `spatalCoords` are not complementary")
  }
  dfPlot <- reshape2::melt(
    as.data.frame(cbind(st.coor, st.pred)), 
    id.vars = c("Spatial 1", "Spatial 2"), 
    variable.name = "CellType", value.name = "Proportion"
  )
  # TODO: change this to parse a gradient palette blablabla
  if (missing(colors)) colors <- scale_color_gradient(low = "white", high = "blue")
  if (is.null(title)) title.plot <- "Predicted proportions"
  plot <- ggplot(
    dfPlot, aes(x = .data[["Spatial 1"]], y = .data[["Spatial 2"]], color = Proportion)
  ) + geom_point(size = size.point, alpha = alpha.point) + colors + 
    ggtitle(title.plot) + 
    SpatialDDLSTheme() 
  if (color.scale == "shared") {
    plot <- plot + facet_wrap(~ CellType, nrow = nrow, ncol = ncol) 
  } 
  return(plot)
}


################################################################################
##################### Plot spatial proportions (single) ########################
################################################################################

#' Generate correlation plots between predicted and expected cell type
#' proportions from test data
#'
#' Generate correlation plot between predicted and expected cell type
#' proportions from test data. Correlation plots can be displayed all mixed or
#' split by cell type (\code{CellType}) or number of cell types present in the
#' samples (\code{nCellTypes}). See the \code{facet.by} argument and examples
#' for more information. Moreover, a user-selected correlation value is
#' displayed as an annotation on the plots. See the \code{corr} argument for
#' details.
#'
#' @param object \code{\linkS4class{SpatialDDLS}} object with
#'   \code{trained.model} slot containing metrics in the
#'   \code{test.deconv.metrics} slot of a
#'   \code{\linkS4class{SpatialDDLSDNN}} object.
#' @param colors Vector of colors to be used. Only vectors with a number of
#'   colors equal to or greater than the levels of \code{color.by} will be
#'   accepted. By default, a custom color list is used.
#' @param facet.by Variable used to display data in different panels. If
#'   \code{NULL}, the plot is not split into different panels. Options are
#'   \code{nCellTypes} (by number of different cell types) and \code{CellType}
#'   (by cell type).
#' @param color.by Variable used to color data. Options are \code{nCellTypes}
#'   and \code{CellType}.
#' @param corr Correlation value displayed as an annotation on the plot.
#'   Available metrics are Pearson's correlation coefficient (\code{'pearson'})
#'   and concordance correlation coefficient (\code{'ccc'}). The argument can be
#'   \code{'pearson'}, \code{'ccc'} or \code{'both'} (by default).
#' @param filter.sc Boolean indicating whether single-cell profiles are filtered
#'   out and only errors associated with pseudo-bulk samples are displayed
#'   (\code{TRUE} by default).
#' @param pos.x.label X-axis position of correlation annotations (0.95 by
#'   default).
#' @param pos.y.label Y-axis position of correlation annotations (0.1 by
#'   default).
#' @param sep.labels Space separating annotations if \code{corr} is equal to
#'   \code{'both'} (0.15 by default).
#' @param size.point Size of points (0.1 by default).
#' @param alpha.point Alpha of points (0.1 by default).
#' @param nrow Number of rows if \code{facet.by} is different from \code{NULL}.
#' @param ncol Number of columns if \code{facet.by} is other than \code{NULL}.
#' @param title Title of the plot.
#' @param theme \pkg{ggplot2} theme.
#' @param ... Additional arguments for the \link[ggplot2]{facet_wrap} function
#'   from \pkg{ggplot2} if \code{facet.by} is not \code{NULL}.
#'
#' @return A ggplot object with the correlation plots between expected and
#'   actual proportions.
#'
#' @export
#'
#' @seealso \code{\link{calculateEvalMetrics}} \code{\link{distErrorPlot}}
#'   \code{\link{blandAltmanLehPlot}} \code{\link{barErrorPlot}}
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(
#'     counts = matrix(
#'       rpois(30, lambda = 5), nrow = 15, ncol = 20,
#'       dimnames = list(paste0("Gene", seq(15)), paste0("RHC", seq(20)))
#'     )
#'   ),
#'   colData = data.frame(
#'     Cell_ID = paste0("RHC", seq(20)),
#'     Cell_Type = sample(x = paste0("CellType", seq(6)), size = 20,
#'                        replace = TRUE)
#'   ),
#'   rowData = data.frame(
#'     Gene_ID = paste0("Gene", seq(15))
#'   )
#' )
#' DDLS <- loadSCProfiles(
#'   single.cell.data = sce,
#'   cell.ID.column = "Cell_ID",
#'   gene.ID.column = "Gene_ID"
#' )
#' probMatrixValid <- data.frame(
#'   Cell_Type = paste0("CellType", seq(6)),
#'   from = c(1, 1, 1, 15, 15, 30),
#'   to = c(15, 15, 30, 50, 50, 70)
#' )
#' DDLS <- generateBulkCellMatrix(
#'   object = DDLS,
#'   cell.ID.column = "Cell_ID",
#'   cell.type.column = "Cell_Type",
#'   prob.design = probMatrixValid,
#'   num.bulk.samples = 50,
#'   verbose = TRUE
#' )
#' # training of DDLS model
#' tensorflow::tf$compat$v1$disable_eager_execution()
#' DDLS <- trainSpatialDDLSModel(
#'   object = DDLS,
#'   on.the.fly = TRUE,
#'   batch.size = 15,
#'   num.epochs = 5
#' )
#' # evaluation using test data
#' DDLS <- calculateEvalMetrics(
#'   object = DDLS
#' )
#' # correlations by cell type
#' corrExpPredPlot(
#'   object = DDLS,
#'   facet.by = "CellType",
#'   color.by = "CellType",
#'   corr = "both"
#' )
#' # correlations of all samples mixed
#' corrExpPredPlot(
#'   object = DDLS,
#'   facet.by = NULL,
#'   color.by = "CellType",
#'   corr = "ccc",
#'   pos.x.label = 0.2,
#'   alpha.point = 0.3
#' )
#' }
#' 
plotSpatialProp <- function(
    object,
    index.st,
    colors,
    cell.type,
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
  if (!cell.type %in% colnames(st.pred)) stop("`cell.type` must be a valid cell type")
  st.pred <- st.pred[, cell.type, drop = FALSE]
  dfPlot <- as.data.frame(cbind(st.coor, st.pred))
  
  if (!all(rownames(st.coor) == rownames(st.pred))) {
    stop("Matrices in `deconv.spots` and `spatalCoords` are not complementary")
  }
  # TODO: change this to parse a gradient palette blablabla
  if (missing(colors)) colors <- scale_color_gradient(low = "white", high = "blue")
  if (is.null(title)) title.plot <- paste0("Predicted proportions (", cell.type, ")")
  
  plot <- ggplot(
    dfPlot, aes(
      x = .data[["Spatial 1"]], y = .data[["Spatial 2"]], 
      color = .data[[cell.type]]
    )
  ) + geom_point(size = size.point, alpha = alpha.point) + colors + 
    ggtitle(title.plot) + SpatialDDLSTheme() 

  return(plot)
}

################################################################################
##################### Plot spatial proportions (blended) #######################
################################################################################

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
  if (!all(rownames(st.coor) == rownames(st.pred))) {
    stop("Matrices in `deconv.spots` and `spatalCoords` are not complementary")
  }
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
    # ggnewscale::new_scale("fill") +
    ggnewscale::new_scale_color() +
    geom_point(
      aes(fill = .data[[cell.type[2]]]),  color = "black",
      size = size.point, alpha = alpha.point, shape = 21
    ) + 
    # scale_color_gradientn(cell.type[2], colors = c("white", "red")) + 
    scale_fill_gradientn(cell.type[2], colors = c("white", "red")) + 
    SpatialDDLSTheme() 
  
  return(plot)
}
