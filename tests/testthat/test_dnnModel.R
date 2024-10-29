context("Training of deconvolution models (deep neural networks): dnnModel.R")

skip_if_not(.checkPythonDependencies(alert = "none"))

# to make compatible with any computer: disable eager execution
tensorflow::tf$compat$v1$disable_eager_execution()

################################################################################
########################## trainDeconvModel function ###########################
################################################################################

# simulating data
set.seed(123)
sce <- SingleCellExperiment::SingleCellExperiment(
  matrix(
    stats::rpois(100, lambda = 5), nrow = 40, ncol = 30, 
    dimnames = list(paste0("Gene", seq(40)), paste0("RHC", seq(30)))
  ),
  colData = data.frame(
    Cell_ID = paste0("RHC", seq(30)),
    Cell_Type = sample(x = paste0("CellType", seq(4)), size = 30, replace = TRUE)
  ),
  rowData = data.frame(
    Gene_ID = paste0("Gene", seq(40))
  )
)
SDDLS <- createSpatialDDLSobject(
  sc.data = sce,
  sc.cell.ID.column = "Cell_ID",
  sc.gene.ID.column = "Gene_ID",
  sc.filt.genes.cluster = FALSE
)
SDDLS <- genMixedCellProp(
  object = SDDLS,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_Type",
  num.sim.spots = 30,
  verbose = TRUE
)
SDDLS <- estimateZinbwaveParams(
  object = SDDLS,
  cell.type.column = "Cell_Type",
  cell.ID.column = "Cell_ID",
  gene.ID.column = "Gene_ID",
  verbose = FALSE
)
SDDLS <- simSCProfiles(
  object = SDDLS,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_Type",
  n.cells = 15,
  verbose = FALSE
)

# simulating spatial data
simSpatialExperiment <- function(n = 1) {
  sim.samples <- function() {
    ngenes <- sample(5:40, size = 1)
    ncells <- sample(5:40, size = 1)
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

# check if object contains all information needed
test_that(
  "Wrong object: lack of specific data", 
  {
    # object without prob.cell.types slot
    expect_error(
      trainDeconvModel(object = SDDLS, num.epochs = 10, verbose = FALSE), 
      regexp = "If `type.data.train` = mixed is selected, 'mixed.profiles' must be provided"
    )
    SDDLS <- suppressWarnings(
      genMixedCellProp(
        object = SDDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        num.sim.spots = 100,
        verbose = FALSE
      )
    )
    # combine = 'both' without mixed samples
    expect_error(
      trainDeconvModel(
        object = SDDLS, type.data.train = "both", 
        num.epochs = 10, verbose = FALSE
      ), 
      regexp = "If `type.data.train` = both'"
    )
    expect_error(
      trainDeconvModel(
        object = SDDLS, type.data.train = "single-cell", 
        num.epochs = 10, type.data.test = "both", verbose = FALSE
      ), 
      regexp = "If `type.data.test` = both' is selected, 'mixed.profiles' and at least one single cell slot must be provided"
    )
    # combine = 'mixed' without mixed samples
    expect_error(
      trainDeconvModel(
        object = SDDLS, type.data.train = "mixed", 
        num.epochs = 10, verbose = FALSE
      ), 
      regexp = "If `type.data.train` = mixed is selected, 'mixed.profiles' must be provided"
    )
    # combine = 'single-cell' without mixed for test data (evaluation of the model)
    expect_error(
      trainDeconvModel(
        object = SDDLS, type.data.train = "single-cell", verbose = FALSE
      ), 
      regexp = "If `type.data.test` = mixed is selected, 'mixed.profiles' must be provided"
    )
    # combine = 'single-cell' without mixed data for test data when on.the.fly = TRUE
    expect_message(
      suppressWarnings(
        trainDeconvModel(
          object = SDDLS,
          type.data.train = "single-cell",
          on.the.fly = TRUE,
          batch.size = 12,
          num.epochs = 10, 
          view.metrics.plot = FALSE
        )
      ),
      regexp = "Training and test on the fly was selected"
    )
    SDDLSBad <- SDDLS
    mixed.profiles(SDDLSBad) <- NULL
    trained.model(SDDLSBad) <- NULL
    # type.data.train = 'mixed' without mixed for test data when on.the.fly = TRUE
    expect_message(
      SDDLSBad <- suppressWarnings(
        trainDeconvModel(
          object = SDDLSBad,
          type.data.train = "mixed",
          on.the.fly = TRUE,
          num.epochs = 10, 
          batch.size = 12,
          view.metrics.plot = FALSE
        )
      ), 
      regexp = "Training and test on the fly was selected"
    )
    # check function to generate pseudo-bulk samples
    expect_error(
      trainDeconvModel(
        object = SDDLS,
        on.the.fly = TRUE,
        agg.function = "Invalid",
        num.epochs = 10, 
        batch.size = 12,
        view.metrics.plot = FALSE
      ), 
      regexp = "'agg.function' must be one of the following options"
    )
  }
)

# check expected parameters
test_that(
  desc = "Parameters", 
  code = {
      SDDLS <- suppressWarnings(
        genMixedCellProp(
          object = SDDLS,
          cell.ID.column = "Cell_ID",
          cell.type.column = "Cell_Type",
          num.sim.spots = 100,
          verbose = FALSE
        )
      )
      SDDLS <- simMixedProfiles(SDDLS, verbose = FALSE)
      # change neural network architecture
      SDDLS <- trainDeconvModel(
        object = SDDLS,
        num.hidden.layers = 3,
        num.units = c(200, 200, 100),
        num.epochs = 10, 
        batch.size = 20,
        verbose = FALSE
      )
      expect_true(
        grepl(
          pattern = "Dense3", 
          as.character(keras::get_config(trained.model(SDDLS)@model))
        )
      )
      trained.model(SDDLS) <- NULL
      SDDLS <- trainDeconvModel(
        object = SDDLS,
        num.hidden.layers = 1,
        num.units = c(100),
        batch.size = 20,
        num.epochs = 10, 
        verbose = FALSE
      )
      expect_false(
        grepl(
          pattern = "Dense3", 
          as.character(keras::get_config(trained.model(SDDLS)@model))
        )
      )
      expect_false(
        grepl("200", as.character(keras::get_config(trained.model(SDDLS)@model)))
      )
      expect_true(
        grepl("100", as.character(keras::get_config(trained.model(SDDLS)@model)))
      )
      # incorrect architecture
      trained.model(SDDLS) <- NULL
      expect_error(
        trainDeconvModel(
          object = SDDLS,
          num.hidden.layers = 1,
          num.units = c(200, 200, 100),
          batch.size = 20,
          num.epochs = 10, 
          verbose = FALSE
        ),
        regexp = "The number of hidden layers must be equal"
      )
      # check if activation.fun works
      SDDLS <- trainDeconvModel(
        object = SDDLS,
        num.hidden.layers = 1,
        num.units = c(100),
        activation.fun = "elu",
        batch.size = 20,
        num.epochs = 10, 
        verbose = FALSE
      )
      expect_true(
        grepl("elu", as.character(keras::get_config(trained.model(SDDLS)@model)))
      )
      expect_false(
        grepl("relu", as.character(keras::get_config(trained.model(SDDLS)@model)))
      )
      # check if dropout.rate works
      trained.model(SDDLS) <- NULL
      SDDLS <- trainDeconvModel(
        object = SDDLS,
        num.hidden.layers = 2,
        num.units = c(100, 100),
        dropout.rate = 0.45,
        batch.size = 20,
        num.epochs = 10, 
        verbose = FALSE
      )
      expect_true(
        grepl("0.45", as.character(keras::get_config(trained.model(SDDLS)@model)))
      )
      # check if loss and metrics work
      trained.model(SDDLS) <- NULL
      SDDLS <- trainDeconvModel(
        object = SDDLS,
        num.hidden.layers = 2,
        num.units = c(100, 100),
        loss = "mean_squared_error",
        metrics = c("accuracy", "mean_absolute_error",
                    "cosine_similarity"),
        batch.size = 20,
        num.epochs = 10, 
        verbose = FALSE
      )
      expect_true(
        any(grepl("accuracy", names(trained.model(SDDLS)@test.metrics)))
      )
      expect_true(
        any(grepl("cosine_similarity", names(trained.model(SDDLS)@test.metrics)))
      )
      # check scaling parameters
      expect_error(
        object = trainDeconvModel(
          object = SDDLS,
          batch.size = 20,
          scaling = "invalid",
          verbose = FALSE
        ), 
        regexp = "'scaling' argument must be one of the following options"
      )
      # check behaviour
      SDDLS@trained.model <- NULL
      SDDLS.standardize <- trainDeconvModel(
        object = SDDLS,
        batch.size = 20,
        scaling = "standardize",
        num.epochs = 10, 
        verbose = FALSE
      )
      SDDLS.rescale <- trainDeconvModel(
        object = SDDLS,
        batch.size = 20,
        scaling = "rescale",
        num.epochs = 10, 
        verbose = FALSE
      )
      samp.stand <- as.numeric(
        gsub(
          pattern = "test_|Spot_|CellType\\d\\_Simul|RH|C", 
          replacement = "", 
          x = rownames(SDDLS.standardize@trained.model@test.pred)
        )
      ) %>% order()
      samp.rescale <- as.numeric(gsub(
        pattern = "test_|Spot_|CellType\\d\\_Simul|RH|C", 
        replacement = "", 
        x = rownames(SDDLS.rescale@trained.model@test.pred)
      )) %>% order()
      stand <- SDDLS.standardize@trained.model@test.pred[samp.stand, ]
      rescale <- SDDLS.rescale@trained.model@test.pred[samp.rescale, ]
      expect_false(object = all(stand == rescale))
    }
)

# check custom.model parameter
test_that(
  desc = "custom.model parameter", 
  {
    SDDLS <- suppressWarnings(
      genMixedCellProp(
        object = SDDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        num.sim.spots = 100,
        verbose = FALSE
      )
    )
    SDDLS <- suppressWarnings(simMixedProfiles(SDDLS, verbose = FALSE))
    # 2 hidden layers without dropouts
    customModel <- keras_model_sequential(name = "CustomModel") %>% 
      layer_dense(
        units = 250, 
        input_shape = nrow(single.cell.real(SDDLS)),
        name = "DenseCustom1"
      ) %>% 
      layer_batch_normalization(name = "CustomBatchNormalization1") %>%
      layer_activation(activation = "elu", name = "ActivationELu1") %>% 
      layer_dense(
        units = 150, 
        name = "DenseCustom2"
      ) %>% 
      layer_batch_normalization(name = "CustomBatchNormalization2") %>%
      layer_activation(activation = "elu", name = "ActivationELu2") %>% 
      layer_dense(
        units = ncol(prob.cell.types(SDDLS, "train") %>% prob.matrix()),
        name = "Dense3"
      ) %>% layer_batch_normalization(name = "CustomBatchNormalization3") %>% 
      layer_activation(activation = "softmax", name = "ActivationSoftmax")
    # check is everything works
    SDDLS <- trainDeconvModel(
      object = SDDLS, 
      custom.model = customModel,
      batch.size = 20,
      num.epochs = 10, 
      verbose = FALSE
    )
    expect_s4_class(
      object = SDDLS, 
      class = "SpatialDDLS"
    )
    expect_true(
      grepl(
        pattern = "CustomBatchNormalization2", 
        as.character(keras::get_config(trained.model(SDDLS)@model))
      )
    )
    # incorrect output units (number of cell types) in custom.model
    customModel <- keras_model_sequential(name = "CustomModel") %>% 
      layer_dense(
        units = 250, 
        input_shape = nrow(single.cell.real(SDDLS)),
        name = "DenseCustom1"
      ) %>% 
      layer_batch_normalization(name = "CustomBatchNormalization1") %>%
      layer_activation(activation = "elu", name = "ActivationELu1") %>% 
      layer_dense(
        units = 150, 
        name = "DenseCustom2"
      ) %>% 
      layer_batch_normalization(name = "CustomBatchNormalization2") %>%
      layer_activation(activation = "elu", name = "ActivationELu2") %>% 
      layer_dense(
        units = 2,
        name = "Dense3"
      ) %>% layer_batch_normalization(name = "CustomBatchNormalization3") %>% 
      layer_activation(activation = "softmax", name = "ActivationSoftmax")
    trained.model(SDDLS) <- NULL
    expect_error(
      trainDeconvModel(
        object = SDDLS, 
        custom.model = customModel,
        batch.size = 20,
        num.epochs = 10, 
        verbose = FALSE
      ), regexp = "The number of neurons of the last layer must be equal"
    )
    # incorrect input units (number of genes) in custom.model
    customModel <- keras_model_sequential(name = "CustomModel") %>% 
      layer_dense(
        units = 250, 
        input_shape = 23,
        name = "DenseCustom1"
      ) %>% 
      layer_batch_normalization(name = "CustomBatchNormalization1") %>%
      layer_activation(activation = "elu", name = "ActivationELu1") %>% 
      layer_dense(
        units = 150, 
        name = "DenseCustom2"
      ) %>% 
      layer_batch_normalization(name = "CustomBatchNormalization2") %>%
      layer_activation(activation = "elu", name = "ActivationELu2") %>% 
      layer_dense(
        units = ncol(prob.cell.types(SDDLS, "train") %>% prob.matrix()),
        name = "Dense3"
      ) %>% layer_batch_normalization(name = "CustomBatchNormalization3") %>% 
      layer_activation(activation = "softmax", name = "ActivationSoftmax")
    trained.model(SDDLS) <- NULL
    expect_error(
      trainDeconvModel(
        object = SDDLS, 
        custom.model = customModel,
        batch.size = 20,
        num.epochs = 10, 
        verbose = FALSE
      ), regexp = "The number of neurons of the first layer must be equal to the number of genes"
    )
    # the last activation function is not softmax
    customModel <- keras_model_sequential(name = "CustomModel") %>% 
      layer_dense(
        units = 250, 
        input_shape = nrow(single.cell.real(SDDLS)),
        name = "DenseCustom1"
      ) %>% 
      layer_batch_normalization(name = "CustomBatchNormalization1") %>%
      layer_activation(activation = "elu", name = "ActivationELu1") %>% 
      layer_dense(
        units = 150, 
        name = "DenseCustom2"
      ) %>% 
      layer_batch_normalization(name = "CustomBatchNormalization2") %>%
      layer_activation(activation = "elu", name = "ActivationELu2") %>% 
      layer_dense(
        units = ncol(prob.cell.types(SDDLS, "train") %>% prob.matrix()),
        name = "Dense3"
      ) %>% layer_batch_normalization(name = "CustomBatchNormalization3") %>% 
      layer_activation(activation = "elu", name = "ActivationElu")
    trained.model(SDDLS) <- NULL
    expect_error(
      trainDeconvModel(
        object = SDDLS, 
        custom.model = customModel,
        batch.size = 20,
        num.epochs = 10, 
        verbose = FALSE
      ), regexp = "In order to get proportions as output, the activation function of the last hidden layer must be 'softmax'"
    )
  }
)

################################################################################
######################### deconvSpatialDDLS function ###########################
################################################################################

# deconvolution of new samples 
test_that(
  "Deconvolution of new samples", 
  {
    SDDLS <- suppressWarnings(
      genMixedCellProp(
        object = SDDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        num.sim.spots = 100,
        verbose = FALSE
      )
    )
    SDDLS <- suppressWarnings(simMixedProfiles(SDDLS, verbose = FALSE))
    # check is everything works
    SDDLS <- trainDeconvModel(
      object = SDDLS,
      batch.size = 20,
      num.epochs = 10, 
      verbose = FALSE
    )
    ste <- simSpatialExperiment(n = 1)
    SDDLS <- loadSTProfiles(
      SDDLS, 
      st.data = ste,
      st.spot.ID.column = "Cell_ID",
      st.gene.ID.column = "Gene_ID",
      st.n.slides = 1
    )
    SDDLS <- deconvSpatialDDLS(object = SDDLS, index.st = 1, pca.space = FALSE)
    # expect_true(names(deconv.spots(SDDLS)) == names(spatial.experiments(SDDLS)))
    expect_true(
      nrow(deconv.spots(SDDLS, 1)[["Regularized"]]) == ncol(spatial.experiments(SDDLS, 1))
    )
    expect_true(
      all(rownames(deconv.spots(SDDLS, 1)[["Regularized"]]) == 
            colnames(spatial.experiments(SDDLS, 1)))
    )
    # index.st does not exist
    expect_error(
      deconvSpatialDDLS(
        object = SDDLS, 
        index.st = "not_exists",
        pca.space = FALSE,
        verbose = FALSE
      ), regexp = "spatial.experiment slot does not contain names, so `index.st` must be an integer"
    )
    # simplify.set: generate a new class from two or more cell types
    deconv.spots(SDDLS) <- NULL
    expect_error(
      deconvSpatialDDLS(
        object = SDDLS,
        index.st = 1,
        pca.space = FALSE,
        simplify.set = list(c("Mc", "M")),
        verbose = FALSE
      ), 
      regexp = "Each element in the list must contain the corresponding new class as name"
    )
    deconv.spots(SDDLS) <- NULL
    SDDLS <- deconvSpatialDDLS(
      object = SDDLS,
      index.st = 1,
      pca.space = FALSE,
      simplify.set = list(CellTypesNew = c("CellType2", "CellType4")),
      verbose = FALSE
    )
    expect_type(deconv.spots(SDDLS, 1), type = "list")
    expect_identical(names(deconv.spots(SDDLS, 1)), c("raw", "simpli.set"))
    expect_true(
      any(colnames(deconv.spots(SDDLS, index.st = 1)[["simpli.set"]]) == 
            "CellTypesNew")
    )
    expect_false(
      any(colnames(deconv.spots(SDDLS, index.st = 1)[["raw"]][["Regularized"]]) == 
            "CellTypesNew")
    )
    deconv.spots(SDDLS) <- NULL
    SDDLS <- deconvSpatialDDLS(
      object = SDDLS,
      index.st = 1,
      pca.space = FALSE,
      simplify.set = list(
        CellTypesNew = c("CellType2", "CellType4"), 
        CellTypesNew2 = c("CellType3", "CellType1")
      ), 
      verbose = FALSE
    )
    expect_true(ncol(deconv.spots(SDDLS, index.st = 1)$simpli.set) == 2)
    expect_true(all(
      c("CellTypesNew", "CellTypesNew2") %in% 
        colnames(deconv.spots(SDDLS, index.st = 1)$simpli.set)
    ))
    # simplify.majority: add up proportions to the most abundant cell type
    deconv.spots(SDDLS) <- NULL
    SDDLS <- deconvSpatialDDLS(
      object = SDDLS,
      index.st = 1,
      pca.space = FALSE,
      simplify.majority = list(c("CellType2", "CellType4"), 
                               c("CellType3", "CellType1")),
      verbose = FALSE
    )
    expect_true(
      all(colnames(deconv.spots(SDDLS, index.st = 1)$simpli.maj) == 
            colnames(deconv.spots(SDDLS, index.st = 1)[["raw"]][["Regularized"]]))
    )
    expect_true(
      all(
        names(which(apply(
          X = deconv.spots(SDDLS, index.st = 1)$simpli.maj != 
            deconv.spots(SDDLS, index.st = 1)[["raw"]][["Regularized"]],
          MARGIN = 2,
          FUN = sum) > 0)
        ) %in% c("CellType2", "CellType4", "CellType3", "CellType1")
      )
    )
    # check if both types of simplify can be stored
    deconv.spots(SDDLS) <- NULL
    SDDLS <- deconvSpatialDDLS(
      object = SDDLS,
      index.st = 1,
      pca.space = FALSE,
      simplify.majority = list(c("CellType2", "CellType4"), 
                               c("CellType3", "CellType1")),
      simplify.set = list(
        CellTypesNew = c("CellType2", "CellType4"), 
        CellTypesNew2 = c("CellType3", "CellType1")
      ),
      verbose = FALSE
    )
    expect_true(
      all(names(deconv.spots(SDDLS, 1)) %in% 
            c("raw", "simpli.set", "simpli.majority"))
    )
    barPlotCellTypes(
      data = SDDLS, index.st = 1, colors = default.colors(), set = "simpli.majority"
    )
    barPlotCellTypes(
      data = SDDLS, index.st = 1, colors = default.colors(), set = NULL
    )
    
    ## deconvolution of more than 1 SpatialExperiment object
    ste <- simSpatialExperiment(n = 6) %>% setNames(paste0("ST", 1:6))
    SDDLS <- loadSTProfiles(
      SDDLS, 
      st.data = ste,
      st.spot.ID.column = "Cell_ID",
      st.gene.ID.column = "Gene_ID",
      st.n.slides = 1
    )
    ## pca.space == FALSE because sometimes var is only explained by 1 PC in fake data
    SDDLS <- deconvSpatialDDLS(object = SDDLS, pca.space = FALSE) 
    expect_true(
      all(names(deconv.spots(SDDLS)) == names(spatial.experiments(SDDLS)))
    )
    expect_true(
      nrow(deconv.spots(SDDLS, 4)[["Regularized"]]) == ncol(spatial.experiments(SDDLS, 4))
    )
    expect_true(
      all(rownames(deconv.spots(SDDLS, 1)[["Regularized"]]) == 
            colnames(spatial.experiments(SDDLS, 1)))
    )
    # index.st does not exist
    expect_error(
      deconvSpatialDDLS(
        object = SDDLS, 
        pca.space = FALSE,
        index.st = "no",
        verbose = FALSE
      ), regexp = "`index.st` contains elements not present in spatial.experiments slot"
    )
    # simplify.set: generate a new class from two or more cell types
    deconv.spots(SDDLS) <- NULL
    expect_error(
      deconvSpatialDDLS(
        object = SDDLS,
        index.st = 1,
        pca.space = FALSE,
        simplify.set = list(c("Mc", "M")),
        verbose = FALSE
      ), 
      regexp = "Each element in the list must contain the corresponding new class as name"
    )
  }
)

# check if saving trained models as JSON-like character objects works
test_that(
  desc = "deconvSpatialDDLS: deconvolution of new samples (JSON objects from disk)", 
  {
    SDDLS <- suppressWarnings(
      genMixedCellProp(
        object = SDDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        num.sim.spots = 100,
        verbose = FALSE
      )
    )
    SDDLS <- suppressWarnings(simMixedProfiles(SDDLS, verbose = FALSE))
    SDDLS <- suppressWarnings(
      trainDeconvModel(
        object = SDDLS,
        batch.size = 20,
        num.epochs = 10, 
        verbose = FALSE
      )
    )
    # save SDDLS object as RDS object: transform Python object into a JSON-like character object
    fileTMP <- tempfile()
    saveRDS(object = SDDLS, file = fileTMP)
    # read and check it out
    SDDLSNew <- readRDS(file = fileTMP)
    expect_type(trained.model(SDDLSNew)@model, type = "list")
    expect_type(trained.model(SDDLS)@model, type = "closure")
    expect_s3_class(
      trained.model(SDDLS)@model, class = "keras.engine.sequential.Sequential"
    )
    # recompile and use it to deconvolve new samples
    se <- SummarizedExperiment::SummarizedExperiment(
      matrix(
        stats::rpois(100, lambda = sample(seq(4, 10), size = 100, replace = TRUE)), 
        nrow = 40, ncol = 15, 
        dimnames = list(paste0("Gene", seq(40)), paste0("Bulk", seq(15)))
      )
    )
    SDDLS <- loadSTProfiles(
      object = SDDLS, st.data = simSpatialExperiment(n = 1),
      st.spot.ID.column = "Cell_ID",
      st.gene.ID.column = "Gene_ID",
      st.n.slides = 1
    )
    SDDLS <- deconvSpatialDDLS(
      object = SDDLS, index.st = 1, verbose = FALSE, pca.space = FALSE
    )
    SDDLSNew <- loadSTProfiles(
      object = SDDLSNew, 
      st.data = simSpatialExperiment(n = 1),
      st.spot.ID.column = "Cell_ID",
      st.gene.ID.column = "Gene_ID",
      st.n.slides = 1
    )
    SDDLSNew <- deconvSpatialDDLS(
      object = SDDLSNew, index.st = 1, pca.space = FALSE, verbose = FALSE
    )
    expect_true(
      all(colnames(deconv.spots(SDDLSNew, 1)) == 
            colnames(deconv.spots(SDDLS, 1)))
    )
    # save DigitalDLSorterDNN object independently of DigitalDLSorter
    fileTMP <- tempfile()
    trainedModelSDDLS <- trained.model(SDDLS)
    saveRDS(object = trainedModelSDDLS, file = fileTMP)
    trainedModelSDDLSNew <- readRDS(file = fileTMP)
    expect_type(model(trainedModelSDDLSNew), type = "list")
    expect_type(model(trainedModelSDDLS), type = "closure")
    expect_s3_class(
      model(trainedModelSDDLS), class = "keras.engine.sequential.Sequential"
    )
  }
)

################################################################################
########################## barPlotCellTypes function ###########################
################################################################################

# visualization of results using barPlotCellTypes function with DigitalDLSorter objects
test_that(
  desc = "barPlotCellTypes: visualization of results using a DigitalDLSorter object", 
  code = {
    SDDLS <- suppressWarnings(
      genMixedCellProp(
        object = SDDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        num.sim.spots = 100,
        verbose = FALSE
      )
    )
    SDDLS <- suppressWarnings(simMixedProfiles(SDDLS, verbose = FALSE))
    # check is everything works
    SDDLS <- suppressWarnings(
      trainDeconvModel(
        object = SDDLS,
        batch.size = 20,
        num.epochs = 10, 
        verbose = FALSE
      )
    )
    SDDLS <- loadSTProfiles(
      SDDLS, 
      st.data = simSpatialExperiment(n = 1),
      st.spot.ID.column = "Cell_ID",
      st.gene.ID.column = "Gene_ID",
      st.n.slides = 1
    )
    SDDLS <- deconvSpatialDDLS(
      object = SDDLS,
      index.st = 1,
      pca.space = FALSE,
      verbose = FALSE
    )
    # index.st not provided
    expect_message(
      barPlotCellTypes(data = SDDLS), 
      regexp = "'index.st' not provided. Setting index.st <- 1"
    )
    # No results available 
    expect_error(
      barPlotCellTypes(data = SDDLS, set = "no_res", index.st = 1),
      regexp = "No simplified results available"
    )
    # invalid simplify argument 
    SDDLS <- deconvSpatialDDLS(
      object = SDDLS,
      index.st = 1,
      pca.space = FALSE,
      simplify.set = list(
        CellTypesNew = c("CellType2", "CellType4"), 
        CellTypesNew2 = c("CellType3", "CellType1")
      ), 
      simplify.majority = list(c("CellType2", "CellType4"), 
                               c("CellType3", "CellType1"))
    )
    expect_error(
      barPlotCellTypes(data = SDDLS, set = "no_res", index.st = 1),
      regexp = "set argument must be one of the following options: 'simpli.set' or 'simpli.majority'"
    )
    # not enough colors
    expect_error(
      barPlotCellTypes(data = SDDLS, colors = c("blue", "red"), index.st = 1),
      regexp = "Number of provided colors is not enough for the number of cell types"
    )
    # incorrect index.st
    expect_error(
      barPlotCellTypes(data = SDDLS, index.st = "no_res"),
      regexp = "Provided 'index.st' does not exist"
    )
    # object without results
    SDDLSBad <- SDDLS
    deconv.spots(SDDLSBad) <- NULL
    expect_error(
      barPlotCellTypes(data = SDDLSBad),
      regexp = "There are no results in SpatialDDLS object."
    )
    # simplify.set and simplify majority work fine --> gg objects
    expect_s3_class(
      barPlotCellTypes(data = SDDLS, index.st = 1), 
      class = "gg"
    )
    expect_s3_class(
      barPlotCellTypes(data = SDDLS, set = "simpli.set", index.st = 1), 
      class = "gg"
    )
    expect_s3_class(
      barPlotCellTypes(data = SDDLS, set = "simpli.majority", index.st = 1), 
      class = "gg"
    )
  }
)
