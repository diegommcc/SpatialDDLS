context("Functions to plot predicted proportions in spatial coordinates: plotSpatialCoor.R")

skip_if_not(.checkPythonDependencies(alert = "none"))

# to make compatible with any computer: disable eager execution
tensorflow::tf$compat$v1$disable_eager_execution()


# simulating data
set.seed(123)
sce <- SingleCellExperiment(
  assays = list(
    counts = matrix(
      stats::rpois(100, lambda = 5), nrow = 40, ncol = 30, 
      dimnames = list(paste0("Gene", seq(40)), paste0("RHC", seq(30)))
    )
  ),
  colData = data.frame(
    Cell_ID = paste0("RHC", seq(30)),
    Cell_Type = sample(x = paste0("CellType", seq(4)), size = 30, replace = TRUE)
  ),
  rowData = data.frame(
    Gene_ID = paste0("Gene", seq(40))
  )
)
simSpatialExperiment <- function(n = 1) {
  sim.samples <- function() {
    ngenes <- sample(3:40, size = 1)
    ncells <- sample(3:40, size = 1)
    counts <- matrix(
      rpois(ngenes * ncells, lambda = 5), ncol = ncells,
      dimnames = list(paste0("Gene", seq(ngenes)), paste0("Spot", seq(ncells)))
    )
    coordinates <- matrix(
      rep(c(1, 2), ncells), ncol = 2
    )
    return(
      SpatialExperiment::SpatialExperiment(
        assays = list(counts = as.matrix(counts)),
        rowData = data.frame(Gene_ID = paste0("Gene", seq(ngenes))),
        colData = data.frame(Cell_ID = paste0("Spot", seq(ncells))),
        spatialCoords = coordinates
      )
    )
  }
  return(replicate(n = n, expr = sim.samples()))
}
SDDLS <- createSpatialDDLSobject(
  sc.data = sce,
  sc.cell.ID.column = "Cell_ID",
  sc.gene.ID.column = "Gene_ID",
  st.data = simSpatialExperiment(n = 10),
  st.spot.ID.column = "Cell_ID",
  st.gene.ID.column = "Gene_ID"
)
SDDLS <- estimateZinbwaveParams(
  object = SDDLS,
  cell.type.column = "Cell_Type",
  cell.ID.column = "Cell_ID",
  gene.ID.column = "Gene_ID",
  verbose = FALSE
)
# object completed
SDDLSComp <- simSCProfiles(
  object = SDDLS,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_Type",
  n.cells = 15,
  verbose = FALSE
)
SDDLSComp <- genMixedCellProp(
  object = SDDLSComp,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_Type",
  num.sim.spots = 100,
  verbose = FALSE
)
SDDLSComp <- simMixedProfiles(SDDLSComp, verbose = FALSE)
SDDLSComp <- trainDeconvModel(
  object = SDDLSComp,
  batch.size = 28,
  verbose = FALSE
)
SDDLSComp <- calculateEvalMetrics(SDDLSComp)
SDDLSComp <- deconvSpatialDDLS(
  SDDLSComp, simplify.set = list(CellTypesNew = c("CellType2", "CellType4"))
)

# plotSpatialPropAll
test_that(
  desc = "plotSpatialPropAll function", 
  code = {
    expect_error(
      plotSpatialPropAll(object = SDDLS), 
      regexp = "Either predictions"
    )
    expect_error(
      plotSpatialPropAll(object = SDDLSComp, index.st = 20), 
      regexp = "subscript out of bounds"
    )
    # incorrect set
    expect_error(
      plotSpatialPropAll(object = SDDLSComp, index.st = 1, set = "ww"), 
      regexp = "`set`must be one of the following options: 'raw', 'simplify.set', 'simpli.majority'"
    )
  }
)

# plotSpatialProp
test_that(
  desc = "plotSpatialProp function", 
  code = {
    expect_error(
      plotSpatialProp(object = SDDLS), 
      regexp = "Either predictions"
    )
    expect_error(
      plotSpatialProp(object = SDDLSComp, index.st = 20), 
      regexp = "subscript out of bounds"
    )
    # incorrect set
    expect_error(
      plotSpatialProp(object = SDDLSComp, index.st = 1, set = "ww"), 
      regexp = "`set`must be one of the following options: 'raw', 'simplify.set', 'simpli.majority'"
    )
    expect_error(
      plotSpatialProp(object = SDDLSComp, index.st = 1, set = "raw", cell.type = "no"), 
      regexp = "`cell.type` must be a valid cell type"
    )
  }
)
