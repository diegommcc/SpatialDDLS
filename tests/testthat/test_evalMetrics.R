context("Evaluation and metrics: evalMetrics.R")

skip_if_not(.checkPythonDependencies(alert = "none"))

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
SDDLSComp <- suppressWarnings(calculateEvalMetrics(SDDLSComp))

# calculateEvalMetrics
test_that(
  desc = "calculateEvalMetrics function", 
  code = {
    # incorrect object: no trained object
    expect_error(
      calculateEvalMetrics(object = SDDLS), 
      regexp = "The provided object does not have a trained model for evaluation"
    )
    # incorrect object: no prob.cell.types slot
    SDDLSCompBad <- SDDLSComp
    prob.cell.types(SDDLSCompBad) <- NULL
    expect_error(
      calculateEvalMetrics(object = SDDLSCompBad), 
      regexp = "The provided object does not contain actual cell proportions in 'prob.cell.types' slot"
    )
    # check if results are properly stored: only MAE
    SDDLSComp <- calculateEvalMetrics(object = SDDLSComp)
    expect_type(trained.model(SDDLSComp) %>% test.deconv.metrics(), type = "list")
    expect_identical(
      names(trained.model(SDDLSComp) %>% test.deconv.metrics()), 
      c("raw", "allData", "filData")
    )
    expect_true(
      all(lapply(
        trained.model(SDDLSComp) %>% test.deconv.metrics(), names
      )$allData == c("MAE", "MSE"))
    )
    expect_true(
      all(lapply(
        trained.model(SDDLSComp) %>% test.deconv.metrics(), names
      )$filData == c("MAE", "MSE"))
    )
    # aggregated results
    expect_identical(
      names(trained.model(SDDLSComp)@test.deconv.metrics[["allData"]][["MAE"]]),
      c("Sample", "CellType", "pBin", "nCellTypes")
    )
    expect_identical(
      names(trained.model(SDDLSComp)@test.deconv.metrics[["filData"]][["MAE"]]),
      c("Sample", "CellType", "pBin", "nCellTypes")
    )
    
    # both metrics: MAE and MSE
    SDDLSComp <- calculateEvalMetrics(object = SDDLSComp)
    expect_type(trained.model(SDDLSComp) %>% test.deconv.metrics(), type = "list")
    expect_identical(
      names(trained.model(SDDLSComp) %>% test.deconv.metrics()), 
      c("raw", "allData", "filData")
    )
    expect_identical(
      lapply(
        trained.model(SDDLSComp) %>% test.deconv.metrics(), names
      )$allData,  c("MAE", "MSE")
    )
    expect_identical(
      lapply(
        trained.model(SDDLSComp) %>% test.deconv.metrics(), names
      )$filData,  c("MAE", "MSE")
    )
    # aggregated results
    expect_identical(
      lapply(trained.model(SDDLSComp)@test.deconv.metrics[["allData"]], names),
      list(
        MAE = c("Sample", "CellType", "pBin", "nCellTypes"), 
        MSE = c("Sample", "CellType", "pBin", "nCellTypes")
      )
    )
    expect_identical(
      lapply(trained.model(SDDLSComp)@test.deconv.metrics[["filData"]], names),
      list(
        MAE = c("Sample", "CellType", "pBin", "nCellTypes"), 
        MSE = c("Sample", "CellType", "pBin", "nCellTypes")
      )
    )
  }
)

# distErrorPlot
test_that(
  desc = "distErrorPlot function", 
  code = {
    # incorrect object: no evaluation metrics
    expect_error(
      distErrorPlot(object = SDDLS, error = "AbsErr"), 
      regexp = "The provided object does not contain evaluation metrics. Use 'calculateEvalMetrics' function"
    )
    # incorrect error parameter
    expect_error(
      distErrorPlot(object = SDDLSComp, error = "no.metrics"), 
      regexp = "'error' provided is not valid"
    )
    # incorrect number of colors
    expect_error(
      distErrorPlot(
        object = SDDLSComp, error = "AbsErr", colors = c("red", "blue")
      ), 
      regexp = "Number of provided colors not enough"
    )
    # incorrect X variable (x.by parameter)
    expect_error(
      distErrorPlot(
        object = SDDLSComp, error = "AbsErr", x.by = "no.variable"
      ), 
      regexp = "'x.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' and 'pBin'"
    )
    # incorrect facet.by parameter
    expect_error(
      distErrorPlot(
        object = SDDLSComp, error = "AbsErr", facet.by = "no.variable"
      ), 
      regexp = "'facet.by' provided is not valid. Available options are: 'nCellTypes', 'CellType' or NULL"
    )
    # incorrect color.by parameter
    expect_error(
      distErrorPlot(
        object = SDDLSComp, error = "AbsErr", color.by = "no.variable"
      ), 
      regexp = "'color.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' and NULL"
    )
    # incorrect type of plot
    expect_error(
      distErrorPlot(
        object = SDDLSComp, error = "AbsErr", type = "no.type"
      ), 
      regexp = "'type' provided is not valid. The available options are: 'violinplot' and 'boxplot'"
    )
    # filtering of single-cell profiles
    p1 <- distErrorPlot(
      object = SDDLSComp, error = "AbsErr", filter.sc = TRUE
    )
    p2 <- distErrorPlot(
      object = SDDLSComp, error = "AbsErr", filter.sc = FALSE
    )
    expect_true(nrow(p1$data) <= nrow(p2$data))
    expect_true(all(grepl(pattern = "Spot", x = p1$data$Sample)))
    expect_false(all(grepl(pattern = "Spot", x = p2$data$Sample)))
  }
)

# corrExpPredPlot
test_that(
  desc = "corrExpPredPlot function", 
  code = {
    # incorrect object: no evaluation metrics
    expect_error(
      corrExpPredPlot(object = SDDLS), 
      regexp = "The provided object does not have evaluation metrics. Use 'calculateEvalMetrics' function"
    )
    # incorrect number of colors
    expect_error(
      corrExpPredPlot(
        object = SDDLSComp, colors = c("red", "blue")
      ), 
      regexp = "The number of provided colors is not enough"
    )
    # incorrect facet.by parameter
    expect_error(
      corrExpPredPlot(
        object = SDDLSComp, facet.by = "no.variable"
      ), 
      regexp = "'facet.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' or NULL"
    )
    # incorrect color.by parameter
    expect_error(
      corrExpPredPlot(
        object = SDDLSComp, color.by = "no.variable"
      ), 
      regexp = "'color.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' or NULL"
    )
    # incorrect correlation
    expect_error(
      corrExpPredPlot(
        object = SDDLSComp, error = "AbsErr", corr = "no.corr"
      ), 
      regexp = "Argument 'corr' invalid. Only supported 'pearson', 'ccc' and 'both'"
    )
    # filtering of single-cell profiles
    p1 <- corrExpPredPlot(object = SDDLSComp, filter.sc = TRUE)
    p2 <- corrExpPredPlot(object = SDDLSComp, filter.sc = FALSE)
    expect_true(nrow(p1$data) <= nrow(p2$data))
    expect_true(all(grepl(pattern = "Spot", x = p1$data$Sample)))
    expect_false(all(grepl(pattern = "Spot", x = p2$data$Sample)))
  }
)

# blandAltmanLehPlot
test_that(
  desc = "blandAltmanLehPlot function", 
  code = {
    # incorrect object: no evaluation metrics
    expect_error(
      blandAltmanLehPlot(object = SDDLS), 
      regexp = "The provided object does not have evaluation metrics. Use 'calculateEvalMetrics' function"
    )
    # incorrect number of colors
    expect_error(
      blandAltmanLehPlot(
        object = SDDLSComp, colors = c("red", "blue")
      ), 
      regexp = "The number of provided colors is not enough"
    )
    # incorrect facet.by parameter
    expect_error(
      blandAltmanLehPlot(
        object = SDDLSComp, facet.by = "no.variable"
      ), 
      regexp = "'facet.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' or NULL"
    )
    # incorrect color.by parameter
    expect_error(
      blandAltmanLehPlot(
        object = SDDLSComp, color.by = "no.variable"
      ), 
      regexp = "'color.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' or NULL"
    )
    # filtering of single-cell profiles
    p1 <- blandAltmanLehPlot(object = SDDLSComp, filter.sc = TRUE)
    p2 <- blandAltmanLehPlot(object = SDDLSComp, filter.sc = FALSE)
    expect_true(nrow(p1$data) <= nrow(p2$data))
    expect_true(all(grepl(pattern = "Spot", x = p1$data$Sample)))
    expect_false(all(grepl(pattern = "Spot", x = p2$data$Sample)))
  }
)


# barErrorPlot
test_that(
  desc = "barErrorPlot function", 
  code = {
    # incorrect object: no evaluation metrics
    expect_error(
      barErrorPlot(object = SDDLS), 
      regexp = "The provided object does not have evaluation metrics. Use 'calculateEvalMetrics' function"
    )
    # incorrect by parameter
    expect_error(
      barErrorPlot(object = SDDLSComp, by = "no.variable"), 
      regexp = "'by' provided is not valid. The available options are: 'nCellTypes', 'CellType'"
    )
    # incorrect error parameter
    expect_error(
      barErrorPlot(object = SDDLSComp, by = "CellType", error = "no.error"), 
      regexp = "'error' provided is not valid. The available errors are: 'MAE', 'MSE'"
    )
    # incorrect dispersion parameter
    expect_error(
      barErrorPlot(
        object = SDDLSComp, by = "CellType", error = "MSE", dispersion = "no.disp"
      ), 
      regexp = "'dispersion' provided is not valid"
    )
  }
)

