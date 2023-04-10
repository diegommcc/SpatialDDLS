context("Simulation of mixed transcriptional profiles: simMixedSpots.R")

################################################################################
######################## genMixedCellProp function #######################
################################################################################

# simulating data
set.seed(123)
sce <- SingleCellExperiment::SingleCellExperiment(
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
SDDLS <- createSpatialDDLSobject(
  sc.data = sce,
  sc.cell.ID.column = "Cell_ID",
  sc.gene.ID.column = "Gene_ID"
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

# check if object contains all information needed
test_that("Wrong object: Wrong column cell type", {
  # incorrect object
  SDDLSBad <- SDDLS
  single.cell.real(SDDLSBad) <- NULL
  expect_error(
    genMixedCellProp(
      object = SDDLSBad,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      num.sim.spots = 200,
      verbose = TRUE
    ), regexp = "'single.cell.real' slot is empty"
  )
  # incorrect column
  expect_error(
    genMixedCellProp(
      object = SDDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "no_column",
      num.sim.spots = 200,
      verbose = TRUE
    ), regexp = "no_column column is not present in cells.metadata"
  )
})

# check proportion.method argument
test_that(
  desc = "Wrong proportion arguments", 
  code = {
    # incorrect number of elements
    expect_error(
      genMixedCellProp(
        object = SDDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        proportion.method = c(0.6, 0.2, 0.5),
        num.sim.spots = 200,
        verbose = TRUE
      ), 
      regexp = "Proportions in `proportion.method` must add up to 1"
    )
    # not add 100
    expect_error(
      genMixedCellProp(
        object = SDDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        proportion.method = c(0.6, 0.2, 0.1),
        num.sim.spots = 200,
        verbose = TRUE
      ), 
      regexp = "Proportions in `proportion.method` must add up to 1"
    )
    # negative numbers
    expect_error(
      genMixedCellProp(
        object = SDDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        proportion.method = c(0.6, 0.2, -0.2),
        num.sim.spots = 200,
        verbose = TRUE
      ), 
      regexp = "Proportions in `proportion.method cannot be less than zero"
    )
    # checking set one method to 0
    expect_is(
      genMixedCellProp(
        object = SDDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        proportion.method = c(0.6, 0.4, 0),
        num.sim.spots = 200,
        verbose = FALSE
      ), 
      class = "SpatialDDLS"
    )
  }
)
# check n.cells
test_that(
  desc = "Check number samples and cells: argument control and expected output",
  code = {
    ## number of cells per mixed transcriptional profile
    SDDLS1 <- genMixedCellProp(
      object = SDDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      num.sim.spots = 200,
      n.cells = 2,
      verbose = TRUE
    )
    expect_equal(prob.cell.types(SDDLS1, "train") %>% cell.names() %>% ncol(), 2)
    
    # num.sim.spots too low to proportions
    expect_error(
      genMixedCellProp(
        object = SDDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        num.sim.spots = 2,
        verbose = TRUE
      ), 
      regexp = "'num.sim.spots' too low compared to 'train.freq.spots'"
    )
    # dim samples <- 1000 (600 train and 400 test) || n.cells <- 250
    SDDLS.1 <- genMixedCellProp(
      object = SDDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      num.sim.spots = 1000,
      train.freq.spots = 1/2,
      n.cells = 250,
      verbose = FALSE
    )
    # number of bulk samples
    # train matrix
    cell.train.matrix <- prob.cell.types(SDDLS.1, "train") %>% prob.matrix
    expect_equal(nrow(cell.train.matrix), 500)
    expect_true(all(rowSums(cell.train.matrix) == 100))
    # test matrix
    cell.test.matrix <- prob.cell.types(SDDLS.1, "test") %>% prob.matrix
    expect_equal(nrow(cell.test.matrix), 500)
    expect_true(all(rowSums(cell.test.matrix) == 100))
    
    # number of cell types
    # train
    n.cell.train <- prob.cell.types(SDDLS.1, "train") %>% cell.names
    expect_equal(dim(n.cell.train), c(500, 250))
    # test
    n.cell.test <- prob.cell.types(SDDLS.1, "test") %>% cell.names
    expect_equal(dim(n.cell.test), c(500, 250))
    # any shared cell between train and test
    expect_false(any(n.cell.train %in% n.cell.test))
    # with random numbers -------------------------------------------------------
    set.seed(123)
    n.cells <- ceiling(runif(n = 1, min = 100, max = 500))
    num.sim.spots <- ceiling(runif(n = 1, min = 200, max = 5000))
    SDDLS.2 <- genMixedCellProp(
      object = SDDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      num.sim.spots = num.sim.spots,
      n.cells = n.cells,
      verbose = FALSE
    )
    # number of bulk samples
    # total matrix
    cell.train.matrix <- prob.cell.types(SDDLS.2, "train") %>% prob.matrix
    cell.test.matrix <- prob.cell.types(SDDLS.2, "test") %>% prob.matrix
    cell.total.matrix <- rbind(cell.train.matrix, cell.test.matrix)
    expect_equal(nrow(cell.total.matrix), num.sim.spots)
    expect_true(all(rowSums(cell.total.matrix) == 100))
    # number of cell types
    # train
    n.cell.train <- prob.cell.types(SDDLS.2, "train") %>% cell.names
    expect_equal(ncol(n.cell.train), n.cells)
    # test
    n.cell.test <- prob.cell.types(SDDLS.2, "test") %>% cell.names
    expect_equal(ncol(n.cell.test), n.cells)
    # any shared cell between train and test
    expect_false(any(n.cell.train %in% n.cell.test))
    # with random numbers and changing proportions ------------------------------
    n.cells <- ceiling(runif(n = 1, min = 100, max = 500))
    num.sim.spots <- ceiling(runif(n = 1, min = 200, max = 5000))
    SDDLS.3 <- genMixedCellProp(
      object = SDDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      num.sim.spots = num.sim.spots,
      proportion.method = c(0.4, 0.5, 0.1),
      n.cells = n.cells,
      verbose = TRUE
    )
    # number of bulk samples
    # total matrix
    cell.train.matrix <- prob.cell.types(SDDLS.3, "train") %>% prob.matrix
    cell.test.matrix <- prob.cell.types(SDDLS.3, "test") %>% prob.matrix
    cell.total.matrix <- rbind(cell.train.matrix, cell.test.matrix)
    expect_equal(nrow(cell.total.matrix), num.sim.spots)
    expect_true(all(rowSums(cell.total.matrix) == 100))
    # number of cell types
    # train
    n.cell.train <- prob.cell.types(SDDLS.3, "train") %>% cell.names
    expect_equal(ncol(n.cell.train), n.cells)
    # test
    n.cell.test <- prob.cell.types(SDDLS.3, "test") %>% cell.names
    expect_equal(ncol(n.cell.test), n.cells)
    # any shared cell between train and test
    expect_false(any(n.cell.train %in% n.cell.test))
  }
)

################################################################################
########################### simMixedProfiles function ###########################
################################################################################

# check if object contains all information needed
test_that(
  desc = "Wrong object: no single-cell profiles or cell type matrix", 
  code = {
    # no prob matrix
    expect_error(
      simMixedProfiles(
        object = SDDLS,
        type.data = "both"
      ), 
      regexp = "'prob.cell.types' slot is empty"
    )
    SDDLS <- genMixedCellProp(
      object = SDDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      num.sim.spots = 100,
      verbose = FALSE
    )
    SDDLSBad <- SDDLS
    single.cell.real(SDDLSBad) <- NULL
    expect_error(
      simMixedProfiles(
        object = SDDLSBad,
        type.data = "both"
      ), 
      regexp = "There are no real single-cell profiles in SpatialDDLS object"
    )
  }
)

SDDLS <- genMixedCellProp(
  object = SDDLS,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_Type",
  num.sim.spots = 10,
  verbose = FALSE
)
SDDLS <- simMixedProfiles(object = SDDLS, type.data = "both", verbose = FALSE)

# check expected behaviour
test_that(
  desc = "Check expected behaviour", 
  code = {
    # incorrect type.data
    expect_error(
      simMixedProfiles(
        object = SDDLS,
        type.data = "incorrect"
      ), 
      regexp = "'type.data' argument must be one of the following options: 'train', 'test' or 'both'"
    )
    # generate only training data
    SDDLSMod <- SDDLS
    mixed.profiles(SDDLSMod) <- NULL
    SDDLSMod <- simMixedProfiles(
      object = SDDLSMod, type.data = "train", verbose = FALSE
    )
    expect_null(mixed.profiles(SDDLSMod, type.data = "test"))
    # dimensions of resulting matrices
    mixed.profiles(SDDLS) <- NULL
    SDDLS <- simMixedProfiles(object = SDDLS, type.data = "both", verbose = FALSE)
    expect_true(
      dim(mixed.profiles(SDDLS, type.data = "train"))[2] == 
        dim(prob.cell.types(SDDLS, type.data = "train")@prob.matrix)[1]
    )
    # chcek how pseudo-bulks are generated
    SDDLS@mixed.profiles <- NULL
    expect_error(
      object = simMixedProfiles(
        object = SDDLS, type.data = "train", 
        mixing.function = "Invalid", verbose = FALSE
      ), 
      regexp = "'mixing.function' must be one of the following options"
    )
    # check if pseudo-bulks are different
    SDDLS.meanCPM <- simMixedProfiles(
      object = SDDLS, type.data = "train", 
      mixing.function = "MeanCPM", verbose = FALSE
    )
    SDDLS.addCPM <- simMixedProfiles(
      object = SDDLS, type.data = "train", 
      mixing.function = "AddCPM", verbose = FALSE
    )
    SDDLS.addCounts <- simMixedProfiles(
      object = SDDLS, type.data = "train", 
      mixing.function = "AddRawCount", verbose = FALSE
    )
    expect_false(
      all(assay(SDDLS.meanCPM@mixed.profiles$train) == 
            assay(SDDLS.addCPM@mixed.profiles$train))
    )
    expect_false(
      all(assay(SDDLS.addCounts@mixed.profiles$train) == 
            assay(SDDLS.addCPM@mixed.profiles$train))
    )
  }
)

# simMixedProfiles: check parameters related to HDF5 files used as back-end
test_that(
  desc = paste(
    "Using HDF5 files as back-end simMixedProfiles: the following", 
    "tests will write file in temp directory/files. Only available if", 
    "DelayedArray and HDF5Array packages are installed"
  ), code = {
    skip_if_not_installed("DelayedArray")
    skip_if_not_installed("HDF5Array")
    mixed.profiles(SDDLS) <- NULL
    # check if HDF5 file exists and if it is correct
    file <- tempfile()
    expect_message(
      SDDLS <- simMixedProfiles(
        object = SDDLS,
        type.data = "both",
        file.backend = file,
        verbose = TRUE
      ), 
      regexp = "=== Writing data to HDF5 file"
    )
    expect_equal(
      dim(mixed.profiles(SDDLS, "train"))[2], 
      dim(prob.cell.types(SDDLS, type.data = "train")@prob.matrix)[1]
    )
    expect_true(file.exists(file))
    expect_s4_class(
      object = mixed.profiles(SDDLS, "train")@assays@data$counts, 
      class = "HDF5Array"
    )
    # file.backend that already exists
    expect_error(
      object = simMixedProfiles(
        object = SDDLS,
        type.data = "both",
        file.backend = file,
        verbose = TRUE
      ),
      regexp = "'file.backend' already exists. Please provide a correct file path"
    )
    # check if block.processing works
    mixed.profiles(SDDLS) <- NULL
    expect_message(
      SDDLS <- simMixedProfiles(
        object = SDDLS,
        type.data = "both",
        file.backend = tempfile(),
        block.processing = TRUE,
        block.size = 3,
        verbose = TRUE
      ), regexp = "Writing block"
    )
  }
)
