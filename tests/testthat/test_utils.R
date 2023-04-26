context("Utils (helper functions): utils.R")

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
SDDLS <- simSCProfiles(
  object = SDDLS,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_Type",
  n.cells = 15,
  verbose = FALSE
)
SDDLSComp <- genMixedCellProp(
  object = SDDLS,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_Type",
  num.sim.spots = 100,
  verbose = FALSE
)
SDDLSComp <- simMixedProfiles(SDDLSComp, verbose = FALSE)
SDDLS@prob.cell.types <- NULL

# getProbMatrix
test_that(
  desc = "getProbMatrix function", 
  code = {
    # incorrect object: no prob.cell.types slot
    expect_error(
      getProbMatrix(SDDLS, type.data = "train"), 
      regexp = "'prob.cell.types' slot is empty"
    )
    # invalid type.data
    expect_error(
      getProbMatrix(SDDLSComp, type.data = "invalid"), 
      regexp = "'type.data' argument must be 'train' or 'test'"
    )
    # no train data
    SDDLSCompBad <- SDDLSComp
    prob.cell.types(SDDLSCompBad, "train") <- NULL
    expect_error(
      getProbMatrix(SDDLSCompBad, type.data = "train"), 
      regexp = "No train data in 'prob.cell.types' slot"
    )
  }
)

# showProbPlot
test_that(
  desc = "showProbPlot function", 
  code = {
    # incorrect object: no prob.cell.types slot
    expect_error(
      showProbPlot(SDDLS, type.data = "train"), 
      regexp = "'prob.cell.types' slot is empty"
    )
    # invalid type.data
    expect_error(
      showProbPlot(SDDLSComp, type.data = "invalid"), 
      regexp = "'type.data' argument must be 'train' or 'test'"
    )
    # invalid set
    expect_error(
      showProbPlot(SDDLSComp, type.data = "train", set = 7), 
      regexp = "'set' argument must be an integer between 1 and 3"
    )
    # invalid type.plot
    expect_error(
      showProbPlot(
        SDDLSComp, type.data = "train", set = 1, type.plot = "no.type"
      ), 
      regexp = "'type.plot' argument must be one of the next options: 'violinplot', 'boxplot', 'linesplot' or 'ncelltypes'"
    )
    # no train data
    SDDLSCompBad <- SDDLSComp
    prob.cell.types(SDDLSCompBad, "train") <- NULL
    expect_error(
      showProbPlot(SDDLSCompBad, type.data = "train", set = 1), 
      regexp = "ProbMatrixCellTypes object does not have plots"
    )
  }
)

# preparingToSave: this is mainly for RDA files
test_that(
  desc = "preparingToSave function", 
  code = {
    skip_if_not(.checkPythonDependencies(alert = "none"))
    SDDLSComp <- trainDeconvModel(
      object = SDDLSComp,
      batch.size = 28,
      verbose = FALSE
    )
    # incorrect object: no trained.model slot
    expect_error(
      preparingToSave(object = SDDLS), 
      regexp = "Provided object has not a DeconvDLModel object"
    )
    # no train data
    SDDLSCompBad <- SDDLSComp
    trained.model(SDDLSCompBad)@model <- list()
    expect_error(
      preparingToSave(object = SDDLSCompBad), 
      regexp = "Provided object has not a trained DNN model"
    )
  }
)

# to keep this variable
fileTMP <- tempfile()

# saving/reading models from/as HDF5 files
test_that(
  desc = "saveTrainedModelAsH5 and loadTrainedModelFromH5: saving/reading models as HDF5 files", 
  code = {
    skip_if_not(.checkPythonDependencies(alert = "none"))
    SDDLSComp <- trainDeconvModel(
      object = SDDLSComp,
      batch.size = 28,
      verbose = FALSE
    )
    # saving model
    # incorrect object: no trained.model slot
    expect_error(
      saveTrainedModelAsH5(object = SDDLS, file.path = fileTMP), 
      regexp = "'trained.model' slot is empty"
    )
    # no train data
    SDDLSCompBad <- SDDLSComp
    trained.model(SDDLSCompBad)@model <- list()
    expect_error(
      saveTrainedModelAsH5(object = SDDLSCompBad, file.path = fileTMP), 
      regexp = "There is not a model to save on disk"
    )
    # save a DNN model from JSON-like character object
    trained.model(SDDLSComp) <- .saveModelToJSON(trained.model(SDDLSComp))
    expect_warning(
      saveTrainedModelAsH5(object = SDDLSComp, file.path = fileTMP), 
      regexp = "Trained model is not a keras object, but a R list with"
    )
    # overwrite a DNN model from JSON-like character object
    expect_warning(
      expect_message(
        saveTrainedModelAsH5(
          object = SDDLSComp, file.path = fileTMP, overwrite = TRUE
        ), 
        regexp = "file already exists. Since 'overwrite' argument is TRUE, it will be overwritten"
      ), regexp = "Trained model is not a keras object"
    )
    # reading model
    # file does not exists
    expect_error(
      loadTrainedModelFromH5(
        object = SDDLSComp, file.path = "no_existent_path"
      ), 
      regexp = "no_existent_path file does not exist"
    )
    # reset.slot = FALSE
    expect_message(
      SDDLSCompNoRes <- loadTrainedModelFromH5(
        object = SDDLSComp, file.path = fileTMP
      ), 
      regexp = "DeconvDLModel object will be overwritten"
    )
    expect_false(is.null(SDDLSCompNoRes@trained.model@training.history))
    # reset.slot = TRUE
    expect_message(
      SDDLSCompRes <- loadTrainedModelFromH5(
        object = SDDLSComp, file.path = fileTMP, reset.slot = TRUE
      ), 
      regexp = "'trained.model' slot is not empty"
    )
    expect_true(is.null(SDDLSCompRes@trained.model@training.history))
  }
)

# plotTrainingHistory
test_that(
  desc = "plotTrainingHistory", 
  code = {
    skip_if_not(.checkPythonDependencies(alert = "none"))
    SDDLSComp <- trainDeconvModel(
      object = SDDLSComp,
      batch.size = 28,
      verbose = FALSE
    )
    # incorrect object: no trained.model slot
    expect_error(
      plotTrainingHistory(object = SDDLS), 
      regexp = "'trained.model' slot is empty"
    )
    SDDLSCompRes <- loadTrainedModelFromH5(
      object = SDDLSComp, file.path = fileTMP, reset.slot = TRUE
    )
    # no training history
    expect_error(
      plotTrainingHistory(object = SDDLSCompRes), 
      regexp = "There is no training history in provided object"
    )
    # incorrect metrics
    # no training history
    expect_error(
      plotTrainingHistory(object = SDDLSComp, metrics = "invalid_metric"), 
      regexp = "None of the given metrics are in the provided object"
    )
  }
)


print("Check Python dependencies")
# reticulate/tensorflow installation
test_that(
  desc = "reticulate and python/tensorflow checks (if dependencies available)", 
  code = {
    skip_if_not(.checkPythonDependencies(alert = "none"))
    expect_true(.isConda())
    expect_true(.isPython())
    expect_true(.isTensorFlow())
  }
)
