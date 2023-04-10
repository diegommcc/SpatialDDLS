context("Simulation of single-cell RNA-Seq profiles: simSingleCell.R")

################################################################################
####################### estimateZinbwaveParams function ########################
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
    Cell_Type = sample(x = paste0("CellType", seq(4)), size = 30, replace = TRUE),
    Patient = sample(x = c("Healthy", "Disease"), size = 30, replace = TRUE)
  ),
  rowData = data.frame(
    Gene_ID = paste0("Gene", seq(40)),
    gene_length = sample(seq(50, 1000), size = 40, replace = TRUE)
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

# estimateZinbwaveParams: check if the function detects errors in parameters
test_that(
  desc = "Wrong parameters in estimateZinbwaveParams", 
  code = {
    # incorrect object
    SDDLSBad <- SDDLS
    single.cell.real(SDDLSBad) <- NULL
    zinb.params(SDDLSBad) <- NULL
    expect_error(
      estimateZinbwaveParams(
        object = SDDLSBad,
        cell.ID.column = "Cell_ID",
        gene.ID.column = "Gene_ID",
        cell.type.column = "Cell_Type",
        cell.cov.columns = "Patient",
        gene.cov.columns = "gene_length",
        set.type = "All",
        threads = 1,
        verbose = TRUE
      ), 
      regexp = "'single.cell.real' slot is empty"
    )
    # incorrect column
    SDDLSMod <- SDDLS
    zinb.params(SDDLSMod) <- NULL
    expect_error(
      estimateZinbwaveParams(
        object = SDDLSMod,
        cell.ID.column = "Cell_ID",
        gene.ID.column = "no_column",
        cell.type.column = "Cell_Type",
        cell.cov.columns = "Patient",
        gene.cov.columns = "gene_length",
        set.type = "All",
        threads = 1,
        verbose = TRUE
      ), 
      regexp = "no_column column is not present in genes.metadata"
    )
    # variable with less than two levels
    sce <- single.cell.real(SDDLSMod) 
    colData(sce)$Patient <- 1
    single.cell.real(SDDLSMod) <- sce
    expect_error(
      estimateZinbwaveParams(
        object = SDDLSMod,
        cell.ID.column = "Cell_ID",
        gene.ID.column = "external_gene_name",
        cell.type.column = "Cell_Type",
        cell.cov.columns = "Patient",
        gene.cov.columns = "gene_length",
        set.type = "All",
        threads = 1,
        verbose = TRUE
      ), 
      regexp = "Patient must have 2 or more unique elements"
    )
    sce <- single.cell.real(SDDLS) 
    rowData(sce)$gene_length <- 1
    single.cell.real(SDDLSMod) <- sce
    expect_error(
      estimateZinbwaveParams(
        object = SDDLSMod,
        cell.ID.column = "Cell_ID",
        gene.ID.column = "external_gene_name",
        cell.type.column = "Cell_Type",
        cell.cov.columns = "Patient",
        gene.cov.columns = "gene_length",
        set.type = "All",
        threads = 1,
        verbose = TRUE
      ), 
      regexp = "gene_length must have 2 or more unique elements"
    )
    # an object with less than two cell types
    SDDLSMod <- SDDLS
    sce <- single.cell.real(SDDLS) 
    colData(sce)$Cell_Type <- 1
    single.cell.real(SDDLSMod) <- sce
    zinb.params(SDDLSMod) <- NULL
    expect_error(
      estimateZinbwaveParams(
        object = SDDLSMod,
        cell.ID.column = "Cell_ID",
        gene.ID.column = "Gene_ID",
        cell.type.column = "Cell_Type",
        cell.cov.columns = "Patient",
        gene.cov.columns = "gene_length",
        set.type = "All",
        threads = 1,
        verbose = TRUE
      ), 
      regexp = "'cell.type.column' must have 2 or more unique elements"
    )
    # incorrect set.type
    SDDLSMod <- SDDLS
    zinb.params(SDDLSMod) <- NULL
    expect_error(
      estimateZinbwaveParams(
        object = SDDLSMod,
        cell.ID.column = "Cell_ID",
        gene.ID.column = "Gene_ID",
        cell.type.column = "Cell_Type",
        cell.cov.columns = "Patient",
        gene.cov.columns = "gene_length",
        set.type = "non-existent-cell-type",
        verbose = FALSE
      ), 
      regexp = "provided in 'set.type' argument not found"
    )
  }
)
# estimateZinbwaveParams: .reduceDataset function
test_that(
  desc = "Functions to subset data in estimateZinwaveParams function (.reduceDataset)", 
  code = {
    list.data <- .extractDataFromSE(
      SEobject = single.cell.real(SDDLS),
      cell.ID.column = "Cell_ID",
      gene.ID.column = "Gene_ID",
      new.data = FALSE
    )
    # incorrect object
    expect_error(
      .reduceDataset(
        subset.cells = 100,
        list.data = list.data,
        cell.type.column = "Cell_Type", 
        cell.ID.column = "Cell_ID", 
        gene.ID.column = "Gene_ID",
        proportional = FALSE,
        verbose = FALSE
      ), 
      regexp = "'subset.cells' must be less than the total number of cells"
    )
    expect_error(
      .reduceDataset(
        subset.cells = 2,
        list.data = list.data,
        cell.type.column = "Cell_Type", 
        cell.ID.column = "Cell_ID", 
        gene.ID.column = "Gene_ID",
        proportional = FALSE,
        verbose = FALSE
      ), 
      regexp = "'subset.cells' must be greater than the number of cell types"
    )
  }
)
# check object ZinbWave
test_that(
  desc = "Check object ZinbWave", 
  code = {
    zinb.params(SDDLS) <- NULL
    SDDLS <- estimateZinbwaveParams(
      object = SDDLS,
      cell.ID.column = "Cell_ID",
      gene.ID.column = "Gene_ID",
      cell.type.column = "Cell_Type",
      cell.cov.columns = "Patient",
      gene.cov.columns = "gene_length",
      set.type = "All",
      threads = 1,
      verbose = TRUE
    )
    expect_s4_class(object = zinb.params(SDDLS), class = "ZinbParametersModel")
  }
)

################################################################################
########################### simSCProfiles function #############################
################################################################################

# check if the function detects errors in parameters
test_that(
  desc = "Wrong parameters in simSCProfiles", 
  code = {
    # incorrect object
    SDDLSBad <- SDDLS
    zinb.params(SDDLSBad) <- NULL
    single.cell.simul(SDDLS) <- NULL
    expect_error(
      simSCProfiles(
        object = SDDLSBad,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        n.cells = 10,
        suffix.names = "_Simul",
        verbose = TRUE
      ), 
      regexp = "'zinb.params' slot is empty."
    )
    # incorrect column
    expect_error(
      simSCProfiles(
        object = SDDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "non_existent_column",
        n.cells = 10,
        suffix.names = "_Simul",
        verbose = TRUE
      ), 
      regexp = "non_existent_column column is not present in cells.metadata"
    )
    # n.cells
    expect_error(
      simSCProfiles(
        object = SDDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        n.cells = 0,
        suffix.names = "_Simul",
        verbose = TRUE
      ), 
      regexp = "'n.cells' must be greater than 0"
    )
    # cell.types 
    expect_error(
      simSCProfiles(
        object = SDDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        n.cells = 10,
        cell.types = "non_existent_column",
        suffix.names = "_Simul",
        verbose = TRUE
      ), 
      regexp = "provided in 'cell.types' not found in ZINB-WaVE model"
    )
  }
)
# simSCProfiles: check if parameters work as expected
test_that(
  desc = "Parameters working as expected in simSCProfiles", 
  code = {
    # suffix.names 
    single.cell.simul(SDDLS) <- NULL
    SDDLS <- simSCProfiles(
      object = SDDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      n.cells = 10,
      suffix.names = "_Suffix",
      verbose = FALSE
    )
    expect_true(
      all(grepl(pattern = "_Suffix", colnames(single.cell.simul(SDDLS))))
    )
    expect_true(any(colnames(colData(single.cell.simul(SDDLS))) == "suffix"))
    # warning if suffix column in cells metadata is going to be overwritten
    colData(single.cell.real(SDDLS))$suffix <- 1
    single.cell.simul(SDDLS) <- NULL
    expect_warning(
      simSCProfiles(
        object = SDDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        n.cells = 10,
        suffix.names = "_Simul",
        verbose = FALSE
      )
    )
    # correct number of cells
    single.cell.simul(SDDLS) <- NULL
    colData(single.cell.real(SDDLS))$suffix <- NULL
    SDDLS <- simSCProfiles(
      object = SDDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      n.cells = 14,
      suffix.names = "_Suffix",
      verbose = FALSE
    )
    expect_equal(dim(single.cell.simul(SDDLS))[2], 14 * 4)
    single.cell.simul(SDDLS) <- NULL
    SDDLS <- simSCProfiles(
      object = SDDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      n.cells = 14,
      suffix.names = "_Suffix",
      cell.types = c("CellType1", "CellType3"),
      verbose = FALSE
    )
    expect_equal(dim(single.cell.simul(SDDLS))[2], 14 * 2)
    expect_true(
      all(
        unique(colData(single.cell.simul(SDDLS))[["Cell_Type"]]) %in% 
          c("CellType1", "CellType3")
      )
    )
  }
)

# simSCProfiles: check parameters related to HDF5 files used as back-end
test_that(
  desc = paste("Using HDF5 files as back-end simSCProfiles: the following", 
        "tests will write file in temp directory/files. Only available if", 
        "DelayedArray and HDF5Array packages are installed"), 
  code = {
    skip_if_not_installed("DelayedArray")
    skip_if_not_installed("HDF5Array")
    # check if HDF5 file exists and if it is correct
    single.cell.simul(SDDLS) <- NULL
    file <- tempfile()
    expect_message(
      SDDLS <- simSCProfiles(
        object = SDDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        n.cells = 10,
        file.backend = file,
        verbose = TRUE
      ), 
      regexp = "=== Writing data to HDF5 file"
    )
    expect_equal(dim(single.cell.simul(SDDLS))[2], 10 * 4)
    expect_true(file.exists(file))
    expect_s4_class(object = counts(single.cell.simul(SDDLS)), class = "HDF5Array")
    # check if name.dataset.backend changes the name of dataset used
    single.cell.simul(SDDLS) <- NULL
    SDDLS <- simSCProfiles(
      object = SDDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      n.cells = 10,
      file.backend = file,
      name.dataset.backend = "new.dataset",
      verbose = FALSE
    )
    expect_true("new.dataset" %in% rhdf5::h5ls(file)[, "name"])
    # cannot be used the same dataset in the same HDF5 file
    single.cell.simul(SDDLS) <- NULL
    expect_error(
      simSCProfiles(
        object = SDDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        n.cells = 10,
        file.backend = file,
        name.dataset.backend = "new.dataset",
        verbose = FALSE
      )
    )
    # check if block.processing works
    single.cell.simul(SDDLS) <- NULL
    expect_message(
      SDDLS <- simSCProfiles(
        object = SDDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        n.cells = 10,
        file.backend = file,
        name.dataset.backend = "new.dataset.1",
        block.processing = TRUE,
        block.size = 5,
        verbose = TRUE
      ), regexp = "=== Simulating and writing new single-cell profiles by block"
    )
  }
)
