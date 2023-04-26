context("Loading scRNA-seq data into the SpatialDDLS object: loadData.R")

################################################################################
##################### From a SingleCellExperiment object #######################
################################################################################

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
# errors related to wrong columns metadata
test_that(
  desc = "Wrong metadata columns return errors", 
  code = {
    expect_error(
      createSpatialDDLSobject(
        sc.data = sce,
        sc.cell.ID.column = "CellID",
        sc.gene.ID.column = 1
      )
    )
    expect_error(
      suppressWarnings(createSpatialDDLSobject(
        sc.sc.data = sce,
        sc.cell.ID.column = "Cell_type",
        sc.gene.ID.column = 2
      ))
    )
    expect_error(
      createSpatialDDLSobject(
        sc.data = sce,
        sc.cell.ID.column = "non_existent_column",
        sc.gene.ID.column = 2
      )
    )
    expect_error(
      createSpatialDDLSobject(
        sc.data = sce,
        sc.cell.ID.column = "Cell_type",
        sc.gene.ID.column = "non_existent_column"
      )
    )
  }
)

# errors related to remove cells or genes (sc.min.cells and sc.min.counts)
test_that(
  desc = "Catch errors related to sc.min.counts and sc.min.cells", 
  code = {
    expect_error(
      createSpatialDDLSobject(
        sc.data = sce,
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 1,
        sc.min.counts = -1,
        sc.min.cells = 1
      ), 
      regexp = "'min.counts' and 'min.cells' must be greater than or equal to zero"
    )
    expect_error(
      createSpatialDDLSobject(
        sc.data = sce,
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 1,
        sc.min.counts = 1,
        sc.min.cells = -1
      ), 
      regexp = "'min.counts' and 'min.cells' must be greater than or equal to zero"
    )
    expect_error(
      createSpatialDDLSobject(
        sc.data = sce,
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 1,
        sc.min.counts = 10e6,
        sc.min.cells = 10
      ),
      regexp = "Resulting count matrix after filtering"
    )
    expect_error(
      createSpatialDDLSobject(
        sc.data = sce,
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 1,
        sc.min.counts = 30,
        sc.min.cells = 440
      ),
      regexp = "Resulting count matrix after filtering"
    )
  }
)

test_that(
  desc = "Check if filtering works as expected", 
  code = {
    counts.real <- assay(sce)
    sc.min.counts <- 6
    sc.min.cells <- 12
    counts <- counts.real[Matrix::rowSums(counts.real > sc.min.counts) >= sc.min.cells, ]
    SDDLSFiltered <- createSpatialDDLSobject(
      sc.data = sce,
      sc.cell.ID.column = "Cell_ID",
      sc.gene.ID.column = 1,
      sc.min.counts = sc.min.counts,
      sc.min.cells = sc.min.cells
    )
    expect_equal(dim(counts), dim(single.cell.real(SDDLSFiltered)))
    expect_equal(
      as.matrix(counts), as.matrix(assay(single.cell.real(SDDLSFiltered)))
    )
  }
)

test_that(
  desc = "Check if counts matrix is a sparse matrix object", 
  code = {
    DDLS1 <- createSpatialDDLSobject(
      sc.data = sce,
      sc.cell.ID.column = "Cell_ID",
      sc.gene.ID.column = 1,
      sc.min.counts = 0,
      sc.min.cells = 0
    )
    DDLS2 <- createSpatialDDLSobject(
      sc.data = sce,
      sc.cell.ID.column = "Cell_ID",
      sc.gene.ID.column = 1,
      sc.min.counts = 6,
      sc.min.cells = 12
    )
    expect_is(assay(single.cell.real(DDLS1)), class = "dgCMatrix")
    expect_is(assay(single.cell.real(DDLS2)), class = "dgCMatrix")
  }
)

# errors related to SingleCellExperiment: data not provided in some slot,
# data provided in other slots, count matrix does not have row and column
# names
test_that(
  desc = "Wrong SingleCellExperiment object", 
  code = {
    counts <- assay(sce)
    # 1 - no rownames neither rowData: genes
    countsNoGenes <- assay(sce)
    rownames(countsNoGenes) <- NULL
    sceLiNoGenes <- SingleCellExperiment(
      assay = list(counts = countsNoGenes),
      colData = colData(sce)
    )
    expect_error(
      createSpatialDDLSobject(
        sc.data = sceLiNoGenes,
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 1,
        sc.min.counts = 0,
        sc.min.cells = 0
      ), 
      regexp = "Count matrix must have rownames"
    )
    # 2 - no colnames neither colData: cells
    countsNoCells <- assay(sce)
    colnames(countsNoCells) <- NULL
    sceLiNoCells <- SingleCellExperiment(
      assay = list(counts = countsNoCells),
      rowData = rowData(sce)
    )
    expect_error(
      createSpatialDDLSobject(
        sc.data = sceLiNoCells,
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 1,
        sc.min.counts = 0,
        sc.min.cells = 0
      ), 
      regexp = "No data provided in colData slot"
    )
    # 3 - no rowData: genes
    sceLiRNoRowData <- SingleCellExperiment(
      assay = list(counts =  assay(sce)),
      colData = colData(sce)
    )
    expect_error(
      createSpatialDDLSobject(
        sc.data = sceLiRNoRowData,
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 1,
        sc.min.counts = 0,
        sc.min.cells = 0
      ), 
      regexp = "No data provided in rowData slot"
    )
    # 4 - no colnames: cells in matrix
    dfCellsMetadata <- colData(sce)
    rownames(dfCellsMetadata) <- NULL
    sceC <- sce
    colnames(sceC) <- NULL
    sceLiNoColNames <- SingleCellExperiment(
      assay = list(counts = assay(sceC)),
      colData = dfCellsMetadata,
      rowData = rowData(sceC)
    )
    expect_error(
      createSpatialDDLSobject(
        sc.data = sceLiNoColNames,
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 1,
        sc.min.counts = 0,
        sc.min.cells = 0
      ), regexp = "Count matrix must have"
    )
    # 5 - No matrix counts
    sceLiNoCounts <- SingleCellExperiment(
      colData = colData(sce),
      rowData = rowData(sce)
    )
    expect_error(
      createSpatialDDLSobject(
        sc.data = sceLiNoCounts,
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 1,
        sc.min.counts = 0,
        sc.min.cells = 0
      ), 
      regexp = "No count data in SingleCellExperiment object"
    )
    # 6 - More than one assay in SingleCellExperiment: warning, no error
    sceLiMoreThanOne <- SingleCellExperiment(
      assay = list(counts = assay(sce), log = log2(assay(sce)  + 1)),
      colData = colData(sce),
      rowData = rowData(sce)
    )
    expect_warning(
      createSpatialDDLSobject(
        sc.data = sceLiMoreThanOne,
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 1,
        sc.min.counts = 0,
        sc.min.cells = 0
      ), 
      regexp = "There is more than one assay, only the first will be used"
    )
  }
)

test_that(
  desc = "Check if createSpatialDDLSobject works as expected", 
  code = {
    DDLS <- createSpatialDDLSobject(
      sc.data = sce,
      sc.cell.ID.column = "Cell_ID",
      sc.gene.ID.column = 1,
      sc.min.counts = 2,
      sc.min.cells = 2
    )
    expect_is(single.cell.real(DDLS), class = "SingleCellExperiment")
    expect_is(single.cell.simul(DDLS), class = "NULL")
  }
)

################################################################################
##################### From files: tsv and sparse matrices ######################
################################################################################

# core functions are the same for files and SCE objects, so the behavior
# should be the same
file.tests <- "../testdata"
files.tsv <- c("counts.tsv", "cellsMetadata.tsv", "genesMetadata.tsv")
files.tsv.gz <- c("counts.tsv.gz", "cellsMetadata.tsv.gz", "genesMetadata.tsv.gz")
files.sparse <- c("sparse_data/matrix.mtx", "cellsMetadata.tsv", "genesMetadata.tsv")

test_that(
  desc = "Check if loading data from tsv files works as expected", 
  code = {
    expect_message(
      createSpatialDDLSobject(
        sc.data = file.path(file.tests, files.tsv),
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 2,
        sc.min.counts = 2,
        sc.min.cells = 2
      ), "- Filtering features"
    )
    expect_error(
      createSpatialDDLSobject(
        sc.data = file.path(file.tests, files.tsv),
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 1
      )
    )
    expect_error(
      suppressWarnings(
        createSpatialDDLSobject(
          sc.data = file.path(file.tests, files.tsv),
          sc.cell.ID.column = "Cell_type",
          sc.gene.ID.column = 2
        )
      )
    )
    expect_error(
      createSpatialDDLSobject(
        sc.data = file.path(file.tests, files.tsv),
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 2,
        sc.min.counts = 10e6,
        sc.min.cells = 10
      )
    )
  }
)

test_that(
  desc = "Check if loading data from tsv.gz files works as expected", 
  code = {
    expect_message(
      createSpatialDDLSobject(
        sc.data = file.path(file.tests, files.tsv.gz),
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 2,
        sc.min.counts = 2,
        sc.min.cells = 2
      ), "- Filtering features"
    )
    expect_error(
      createSpatialDDLSobject(
        sc.data = file.path(file.tests, files.tsv.gz),
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 1
      )
    )
    expect_error(
      suppressWarnings(
        createSpatialDDLSobject(
          sc.data = file.path(file.tests, files.tsv.gz),
          sc.cell.ID.column = "Cell_type",
          sc.gene.ID.column = 2
        )
      )
    )
    expect_error(
      createSpatialDDLSobject(
        sc.data = file.path(file.tests, files.tsv.gz),
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 2,
        sc.min.counts = 100000,
        sc.min.cells = 10
      )
    )
  }
)

test_that(
  desc = "Check if loading data from sparse files works as expected", 
  code = {
    expect_message(
      createSpatialDDLSobject(
        sc.data = file.path(file.tests, files.sparse),
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 2,
        sc.min.counts = 2,
        sc.min.cells = 2
      ), "- Filtering features"
    )
    expect_error(
      createSpatialDDLSobject(
        sc.data = file.path(file.tests, files.sparse),
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 1
      )
    )
    expect_error(
      suppressWarnings(
        createSpatialDDLSobject(
          sc.data = file.path(file.tests, files.sparse),
          sc.cell.ID.column = "Cell_type",
          sc.gene.ID.column = 2
        )
      )
    )
    expect_error(
      createSpatialDDLSobject(
        sc.data = file.path(file.tests, files.sparse),
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 2,
        sc.min.counts = 100000,
        sc.min.cells = 10
      )
    )
  }
)

test_that(
  desc = "Check if objects from different files are equivalent", 
  code = {
    DDLS.tsv <- createSpatialDDLSobject(
      sc.data = file.path(file.tests, files.tsv),
      sc.cell.ID.column = "Cell_ID",
      sc.gene.ID.column = 2,
      sc.min.counts = 0,
      sc.min.cells = 12
    )
    DDLS.tsv.gz <- createSpatialDDLSobject(
      sc.data = file.path(file.tests, files.tsv.gz),
      sc.cell.ID.column = "Cell_ID",
      sc.gene.ID.column = 2,
      sc.min.counts = 0,
      sc.min.cells = 12
    )
    DDLS.sparse <- createSpatialDDLSobject(
      sc.data = file.path(file.tests, files.sparse),
      sc.cell.ID.column = "Cell_ID",
      sc.gene.ID.column = 2,
      sc.min.counts = 0,
      sc.min.cells = 12
    )
    # DDLS objects
    expect_equal(DDLS.tsv, DDLS.tsv.gz)
    expect_equal(DDLS.tsv, DDLS.sparse)
    # Matrices counts
    expect_equal(
      assay(single.cell.real(DDLS.tsv)),
      assay(single.cell.real(DDLS.tsv.gz))
    )
    expect_equal(
      assay(single.cell.real(DDLS.tsv)),
      assay(single.cell.real(DDLS.sparse))
    )
  }
)

# behaviour of functions with bad built files
test_that(
  desc = "Check if loading data from sparse files works as expected", 
  code = {
    files.tsv.gz.bad.1 <- c(
      "counts.bad.tsv.gz", "cellsMetadata.bad.tsv.gz", "genesMetadata.bad.tsv.gz"
    )
    expect_error(createSpatialDDLSobject(
      sc.data = file.path(file.tests, files.tsv.gz.bad.1),
      sc.cell.ID.column = "Cell_ID",
      sc.gene.ID.column = 2,
      sc.min.counts = 0,
      sc.min.cells = 12
    ))
    files.tsv.gz.bad.2 <- c(
      "counts.bad.tsv.gz", "cellsMetadata.tsv.gz", "genesMetadata.tsv.gz"
    )
    expect_error(createSpatialDDLSobject(
      sc.data = file.path(file.tests, files.tsv.gz.bad.2),
      sc.cell.ID.column = "Cell_ID",
      sc.gene.ID.column = 2,
      sc.min.counts = 0,
      sc.min.cells = 12
    ))
    files.tsv.gz.bad.3 <- c(
      "counts.tsv.gz", "cellsMetadata.bad.tsv.gz", "genesMetadata.tsv.gz"
    )
    expect_error(createSpatialDDLSobject(
      sc.data = file.path(file.tests, files.tsv.gz.bad.3),
      sc.cell.ID.column = "Cell_ID",
      sc.gene.ID.column = 2,
      sc.min.counts = 0,
      sc.min.cells = 12
    ))
    files.tsv.gz.bad.4 <- c(
      "counts.tsv.gz", "cellsMetadata.tsv.gz", "genesMetadata.bad.tsv.gz"
    )
    expect_message(
      object = createSpatialDDLSobject(
        sc.data = file.path(file.tests, files.tsv.gz.bad.4),
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 2,
        sc.min.counts = 0,
        sc.min.cells = 12
      ), 
      regexp = "18008 genes have been discarded from genes metadata"
    )
  }
)

# behaviour of functions with hfd5 files
test_that(
  desc = "Check behaviour with HDF5 files", 
  code = {
    skip_if_not_installed("DelayedArray")
    skip_if_not_installed("HDF5Array")
    file <- tempfile()
    expect_message(
      DDLS.tsv <- createSpatialDDLSobject(
        sc.data = file.path(file.tests, files.tsv),
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 2,
        sc.min.counts = 0,
        sc.min.cells = 12,
        sc.file.backend = file
      ), 
      regexp = "=== Writing data to HDF5 file"
    )
    expect_message(
      DDLS.tsv <- createSpatialDDLSobject(
        sc.data = file.path(file.tests, files.tsv),
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = 2,
        sc.min.counts = 0,
        sc.min.cells = 12,
        sc.file.backend = file,
        sc.name.dataset.backend = "other_dataset",
        sc.block.processing = TRUE
      ),
      regexp = "=== Processing data in HDF5 by blocks"
    )
    expect_true(file.exists(file))
    expect_s4_class(
      object = single.cell.real(DDLS.tsv)@assays@data$counts, 
      class = "HDF5Array"
    )
  }
)

################################################################################
#################### Loading spatial transcriptomics data ######################
################################################################################

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

counts <- matrix(
  rpois(5 * 10, lambda = 5), ncol = 10,
  dimnames = list(paste0("Gene", seq(5)), paste0("Spot", seq(10)))
)
coordinates <- matrix(
  rep(c(1, 2), 10), ncol = 2
)
ste <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = as.matrix(counts)),
  rowData = data.frame(Gene_ID = paste0("Gene", seq(5))),
  colData = data.frame(Cell_ID = paste0("Spot", seq(10))),
  spatialCoords = coordinates
)

# taking shared genes
test_that(
  desc = "Wrong metadata columns return errors", 
  code = {
    expect_error(
      createSpatialDDLSobject(
        sc.data = sce,
        sc.cell.ID.column = "CellID",
        sc.gene.ID.column = 1,
        st.data = sce,
        st.spot.ID.column = "Cell_ID",
        st.gene.ID.column = "Gene_ID"
      ), "`st.data` argument must be a SpatialExperiment"
    )
    # 5 shared genes
    expect_message(
      createSpatialDDLSobject(
        sc.data = sce,
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = "Gene_ID",
        st.data = ste,
        st.spot.ID.column = "Cell_ID",
        st.gene.ID.column = "Gene_ID"
      ), "=== Number of shared genes between single-cell and spatial transcriptomics datasets: 5"
    )
    # deal with list of SpatialExperiment objects
    set.seed(123)
    list.ste <- simSpatialExperiment(n = 5)
    expect_message(
      createSpatialDDLSobject(
        sc.data = sce,
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = "Gene_ID",
        st.data = list.ste,
        st.spot.ID.column = "Cell_ID",
        st.gene.ID.column = "Gene_ID"
      ), "=== 5 SpatialExperiment objects provided"
    )
    ## min number of genes
    set.seed(123)
    list.ste <- simSpatialExperiment(n = 20)
    min.genes <- min(sapply(list.ste, function(x) dim(x)[1]))
    expect_message(
      createSpatialDDLSobject(
        sc.data = sce,
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = "Gene_ID",
        st.data = list.ste,
        st.spot.ID.column = "Cell_ID",
        st.gene.ID.column = "Gene_ID",
        st.n.slides = 20
      ), 
      regexp = paste0("with more genes): ", min.genes)
    )
    # max gene
    set.seed(123)
    list.ste <- simSpatialExperiment(n = 20)
    max.genes <- max(sapply(list.ste, function(x) dim(x)[1]))
    expect_message(
      createSpatialDDLSobject(
        sc.data = sce,
        sc.cell.ID.column = "Cell_ID",
        sc.gene.ID.column = "Gene_ID",
        st.data = list.ste,
        st.spot.ID.column = "Cell_ID",
        st.gene.ID.column = "Gene_ID",
        st.n.slides = 1
      ), 
      paste0("with more genes): ", max.genes)
    )
  }
)

# loading ST data after creating the SpatialDDLS object
test_that(
  desc = "Wrong metadata columns return errors", 
  code = {
    SDDLS <- createSpatialDDLSobject(
      sc.data = sce,
      sc.cell.ID.column = "Cell_ID",
      sc.gene.ID.column = 1
    )
    expect_error(
      loadSTProfiles(
        "fake_obj", 
        st.data = ste,
        st.spot.ID.column = "Cell_ID",
        st.gene.ID.column = "Gene_ID"
      ), "The object provided is not of SpatialDDLS class"
    )
    expect_message(
      loadSTProfiles(
        SDDLS, 
        st.data = ste,
        st.spot.ID.column = "Cell_ID",
        st.gene.ID.column = "Gene_ID"
      ), "=== 1 SpatialExperiment objects provided"
    )
    list.ste <- simSpatialExperiment(n = 10)
    expect_message(
      loadSTProfiles(
        SDDLS, 
        st.data = list.ste,
        st.spot.ID.column = "Cell_ID",
        st.gene.ID.column = "Gene_ID"
      ), "=== 10 SpatialExperiment objects provided"
    )
  }
)
