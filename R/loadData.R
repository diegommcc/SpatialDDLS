#' @importFrom utils read.delim 
#' @importFrom stats setNames
#' @importFrom dplyr desc arrange
#' @importFrom utils head tail
NULL

.readTabFiles <- function(file) {
  if (!file.exists(file)) {
    stop(paste(file, "file provided does not exist"))
  }
  if (grepl(pattern = ".tsv", x = file)) {
    if (grepl(pattern = ".tsv.gz$", x = file)) {
      file.obj <- tryCatch(
        expr = read.delim(file = gzfile(file), sep = "\t",
                          header = T, stringsAsFactors = F),
        error = function(err) {
          stop("The provided file contains duplicated rownames/colnames. ", 
                  "Please, provide a correct count matrix")
        },
        warning = function(err) warning(err)
      )
    } else {
      file.obj <- tryCatch(
        expr = read.delim(file = file, sep = "\t", header = T,
                          stringsAsFactors = F),
        error = function(err) {
          stop("The provided file contains duplicated rownames/colnames. ", 
                  "Please, provide a correct count matrix")
        },
        warning = function(err) warning(err)
      )
    }
  } else if (grepl(pattern = ".rds$", x = file)) {
    file.obj <- readRDS(file = file)
  } else {
    stop("File format is not recognizable. Please, look at allowed formats",
         " in ?loadSCProfiles")
  }
  return(file.obj)
}

.useH5backend <- function(
  counts,
  file.backend,
  compression.level = NULL,
  group = "single.cell",
  chunk.dims = NULL,
  sparse = FALSE,
  verbose = TRUE
) {
  if (!requireNamespace("DelayedArray", quietly = TRUE) || 
      !requireNamespace("HDF5Array", quietly = TRUE)) {
    stop("SpatialDDLS provides the possibility of using HDF5 files as back-end
         when data are too big to be located in RAM. It uses DelayedArray, 
         HDF5Array and rhdf5 to do it. Please install both packages to 
         use this functionality")
  } 
  # if (file.exists(file.backend)) {
  #   if (group %in% rhdf5::h5ls(file.backend)[, "name"]) {
  #     stop("'file.backend' and name group already exist. They cannot exist")  
  #   }
  #   # warning("'file.backend' already exists, but ")
  # }
  if (is.null(compression.level)) {
    compression.level <- HDF5Array::getHDF5DumpCompressionLevel()
  } else {
    if (compression.level < 0 || compression.level > 9) {
      stop("'compression.level' must be an integer between 0 (no ", 
           "compression) and 9 (highest and slowest compression). ")
    }
  }
  if (verbose) message("\n=== Writing data to HDF5 file")
  counts <- DelayedArray::DelayedArray(seed = counts)
  # check correct chunk.dims
  if (is.null(chunk.dims)) {
    chunk.dims <- c(nrow(counts), 1)
  } else {
    if (any(chunk.dims > dim(counts))) {
      warning("'chunk.dims' must be equal to or less than data dimensions. ", 
              "Setting default value", call. = FALSE, immediate. = TRUE)
      chunk.dims <- c(nrow(counts), 1)
    }
  }
  if (sparse) {
    counts <- HDF5Array::writeTENxMatrix(
      x = counts,
      filepath = file.backend,
      group = group,
      level = compression.level,
      verbose = verbose
    )  
  } else {
    counts <- HDF5Array::writeHDF5Array( 
      x = counts,
      filepath = file.backend,
      name = group,
      chunkdim = chunk.dims,
      level = compression.level,
      with.dimnames = TRUE,
      verbose = verbose
    )  
  }
  return(counts)
}

## TODO: this function could be substituted by Seurat::Read10X, but with v5
.readCountsFile <- function(
  counts.file, 
  gene.column = 1,
  name.h5 = NULL,
  file.backend = NULL,
  block.processing = FALSE
) {
  if (grepl(pattern = ".tsv|.rds", x = counts.file, ignore.case = FALSE)) {
    counts <- .readTabFiles(file = counts.file)
  } else if (grepl(pattern = ".mtx$", x = counts.file, ignore.case = FALSE)) {
    if (!file.exists(counts.file))
      stop(paste(counts.file, "file not found"))
    base.dir <- dirname(counts.file)
    if (!file.exists(file.path(base.dir, "genes.tsv"))) 
      stop("No 'genes.tsv' file with mtx file")
    if (!file.exists(file.path(base.dir, "barcodes.tsv"))) 
      stop("No 'barcodes.tsv' file with mtx file")
    counts <- Matrix::readMM(counts.file)
    gene.names <- read.delim(file.path(base.dir, "genes.tsv"), header = F,
                             sep = "\t", stringsAsFactors = F)
    rownames(counts) <- gene.names[, gene.column]
    cell.names <- read.delim(file.path(base.dir, "barcodes.tsv"), header = F,
                             sep = "\t", stringsAsFactors = F)
    colnames(counts) <- cell.names$V1
  } else if (grepl(".h5$|.hdf5$", counts.file, ignore.case = FALSE)) {
    if (is.null(name.h5)) {
      stop("If HDF5 file is provided, the name of dataset used must be given in ",
           "'name.h5' argument") 
    } else if (!is.null(file.backend) && block.processing) {
      # hdf5 file will be used as back-end
      counts <- HDF5Array::HDF5Array(filepath = counts.file, name = name.h5)
    } else if (is.null(file.backend)) {
      # file will be loeaded in memory
      counts <- rhdf5::h5read(file = counts.file, name = name.h5)
    }
  } else {
    stop("File format is not recognizable. Please, look at allowed formats",
         " in ?loadSCProfiles")
  }
  return(counts)
}

.createSCEObject <- function(
  counts, 
  cells.metadata, 
  genes.metadata,
  file.backend,
  name.dataset.backend,
  compression.level,
  chunk.dims,
  block.processing,
  verbose
) {
  # could be a check of counts class -> if (is(counts, "HDF5Array"))
  if (!is.null(file.backend) && 
      !class(counts) %in% c(
        "HDF5Matrix", "HDF5Array", "DelayedArray", "DelayedMatrix"
      )) {
    counts <- .useH5backend(
      counts = counts,
      file.backend = file.backend,
      compression.level = compression.level,
      chunk.dims = chunk.dims,
      group = name.dataset.backend,
      verbose = verbose
    )
  } else if (is.null(file.backend)) {
    counts <- Matrix::Matrix(data = counts, sparse = TRUE)
  }
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = cells.metadata,
    rowData = genes.metadata
  )
  return(sce)
}

.checkColumn <- function(metadata, ID.column, type.metadata, arg) {
  tryCatch(expr = ID.column <- as.numeric(ID.column),
           error = function(e) invisible(x = NULL),
           warning = function(e) invisible(x = NULL))
  if (is(ID.column, "numeric") || is(ID.column, "integer")) {
    if (!ID.column %in% seq(ncol(metadata))) {
      stop(paste(ID.column, "column number is not present in", type.metadata))
    }
  } else if (is(ID.column, "character")) {
    if (!ID.column %in% colnames(metadata)) {
      stop(paste(ID.column, "column is not present in", type.metadata))
    }
  } else {
    stop(paste(arg, "argument is not recognizable"))
  }
}

.processData <- function(
  counts, 
  cells.metadata, 
  cell.ID.column,
  cell.type.column,
  genes.metadata, 
  gene.ID.column,
  filt.genes.cluster,
  min.mean.counts,
  filt.genes.cells,
  min.counts, 
  min.cells,
  file.backend,
  block.processing,
  verbose
) {
  # check if IDs given exist in metadata
  .checkColumn(
    metadata = cells.metadata,
    ID.column = cell.ID.column,
    type.metadata = "cells.metadata",
    arg = "cell.ID.column"
  )
  if (filt.genes.cluster) {
    .checkColumn(
      metadata = cells.metadata,
      ID.column = cell.type.column,
      type.metadata = "cells.metadata",
      arg = "sc.cell.type.column"
    )
  }
  .checkColumn(
    metadata = genes.metadata,
    ID.column = gene.ID.column,
    type.metadata = "genes.metadata",
    arg = "gene.ID.column"
  )
  # duplicated ID cells --------------------------------------------------------
  if (any(duplicated(cells.metadata[, cell.ID.column]))) {
    warning("There are duplicated IDs in 'cells.metadata' (column ", 
            cell.ID.column, "). Making unique")
    cells.metadata[, cell.ID.column] <- make.unique(
      names = cells.metadata[, cell.ID.column]
    )
  }
  # intersect between cells ----------------------------------------------------
  if (!is.null(colnames(counts))) {
    common.cells <- intersect(colnames(counts), cells.metadata[, cell.ID.column])
    diff <- abs(dim(counts)[2] - length(common.cells))
    disc <- abs(length(cells.metadata[, cell.ID.column]) - length(common.cells))
    if (length(common.cells) < min(dim(counts)[2], dim(cells.metadata)[1])) {
      stop(paste("There are ", diff,
                 " cells that don't match between count matrix and metadata"))
    } else if (diff != 0) { # this check includes the previous one
      warning("There are", diff, "cells that don't match between counts ", 
              "matrix and metadata")
    } else if (disc != 0) {
      if (verbose) {
        message("=== Intersection between count matrix and cells metadata:")
        message(
          paste("   ", disc, "cells have been discarded from cells metadata"),
          "\n"
        )
      }
    }
    cells.metadata <- cells.metadata[cells.metadata[, cell.ID.column] %in%
                                       common.cells, , drop = FALSE]
    counts <- counts[, common.cells]
  } else {
    if (ncol(counts) != nrow(cells.metadata)) {
      stop("Count matrix does not have colnames and cells metadata does not ", 
           "have the same number of IDs. Please, provide a correct count matrix")
    } else {
      colnames(counts) <- cells.metadata[, cell.ID.column]
      warning(paste("Count matrix does not have colnames, so", cell.ID.column, 
                    "column of cells metadata will be used")) 
    }
  }
  # intersect between genes ----------------------------------------------------
  if (!is.null(rownames(counts))) {
    common.genes <- intersect(rownames(counts), genes.metadata[, gene.ID.column])
    diff <- abs(dim(counts)[1] - length(common.genes))
    disc <- abs(length(genes.metadata[, gene.ID.column]) - length(common.genes))
    if (length(common.genes) < min(dim(counts)[1], dim(genes.metadata)[1])) {
      stop(paste(
        "There are", diff, 
        "genes that don't match between count matrix and metadata"
      ))
    } else if (diff != 0){
      stop(paste(
        "There are", diff,
        "genes that don't match between count matrix and metadata"
      ))
    } else if (disc != 0) {
      if (verbose) {
        message("    - Intersection between count matrix and genes metadata:")
        message("    ", disc, " genes have been discarded from genes metadata",
                "\n") 
      }
    }
    genes.metadata <- genes.metadata[genes.metadata[, gene.ID.column] %in%
                                       common.genes, , drop = FALSE]
    counts <- counts[common.genes, ]
  } else {
    if (nrow(counts) != nrow(genes.metadata)) {
      stop("Count matrix has not rownames and genes metadata has not the same ", 
           "number of IDs. Please, provide a correct count matrix")
    } else {
      rownames(counts) <- genes.metadata[, gene.ID.column]
      warning(paste("Count matrix has not rownames, so", gene.ID.column, 
                    "column from genes metadata will be used")) 
    } 
  }
  ######### it is here where I can modify the behaviour of the function. Previous 
  ### steps are kind of needed
  
  # filter genes by min.counts and min.cells -----------------------------------
  if (!block.processing) {
    filtered.genes <- .filterGenesSparse(
      counts = counts,
      genes.metadata = genes.metadata,
      gene.ID.column = gene.ID.column,
      cells.metadata = cells.metadata,
      cell.type.column = cell.type.column,
      filt.genes.cluster = filt.genes.cluster,
      min.mean.counts = min.mean.counts,
      filt.genes.cells = filt.genes.cells,
      min.counts = min.counts,
      min.cells = min.cells,
      verbose = verbose
    )  
  } else {
    filtered.genes <- .filterGenesHDF5(
      counts = counts,
      genes.metadata = genes.metadata,
      gene.ID.column = gene.ID.column,
      cells.metadata = cells.metadata,
      cell.type.column = cell.type.column,
      filt.genes.cells = filt.genes.cells,
      min.counts = min.counts,
      min.cells = min.cells,
      filt.genes.cluster = filt.genes.cluster,
      min.mean.counts = min.mean.counts,
      verbose = verbose
    )
  }
  return(list(filtered.genes[[1]], cells.metadata, filtered.genes[[2]]))
}

## solution for large sparse matrices (only works on sparse matrices)
.logicalFiltSparse <- function(counts, min.counts, min.cells) {
  counts <- as(counts, "dgTMatrix") # TsparseMatrix
  dfSumm <- data.frame(
    i = counts@i,
    j = counts@j,
    x = counts@x
  )
  dfSumm <- dfSumm[which(dfSumm$x > min.counts), ]
  ## in case there are no genes that meet the cutoff
  if (any(dim(dfSumm) == 0)) return(rep(FALSE, nrow(counts)))
  m <- Matrix::sparseMatrix(
    i = dfSumm$i + 1,
    j = dfSumm$j + 1,
    x = 1L
  )
  return(Matrix::rowSums(m) >= min.cells)
  # summ <- summary(counts)
  # summ <- summ[which(summ$x > min.counts), ]
  # m <- Matrix::sparseMatrix(
  #   i = summ$i,
  #   j = summ$j,
  #   x = 1L
  # )
}

.filterGenesByCells <- function(
  counts, 
  genes.metadata, 
  min.counts, 
  min.cells
) {
  if (min.counts == 0 && min.cells == 0) {
    return(list(counts, genes.metadata))
  } else if (min.counts < 0 || min.cells < 0) {
    stop("'min.counts' and 'min.cells' must be greater than or equal to zero")
  }
  if (is(counts, "matrix")) {
    counts <- counts[Matrix::rowSums(counts > min.counts) >= min.cells, ]
  } else if (is(counts, "Matrix")) {
    counts <- counts[.logicalFiltSparse(counts, min.counts, min.cells), ]  
  }
  
  return(list(counts, genes.metadata))
}

# TODO: control what numbers can be used for min.mean.counts
# .filterGenesByCluster <- function(
#   counts, 
#   genes.metadata, 
#   cells.metadata, 
#   cell.type.column, 
#   min.mean.counts
# ) {
#   ## mean counts per cluster: with big matrices, this will be problematic
#   # mean.cluster <- aggregate(
#   #   x = t(as.matrix(counts)), 
#   #   by = list(cells.metadata[[cell.type.column]]), 
#   #   FUN = mean
#   # )
#   # rownames(mean.cluster) <- mean.cluster[, 1]
#   # mean.cluster[, 1] <- NULL
#   ## this part has been changed using the same implementation as Matrix.utils
#   mean.cluster <- round(
#     x = .aggregate.Matrix.sparse(
#       x = Matrix::t(counts), 
#       groupings = list(cells.metadata[[cell.type.column]]), 
#       fun = "mean"
#     ), digits = 2
#   ) 
#   ## why this round here??
#   ## using cutoffs
#   sel.genes <- unlist(
#     apply(
#       X = mean.cluster >= min.mean.counts, 
#       MARGIN = 2, 
#       FUN = function(x) if(any(x)) return(TRUE)
#     )
#   )
#   if (!any(sel.genes)) {
#     return(list(counts, genes.metadata))
#   } else {
#     return(
#       list(
#         counts[names(sel.genes), ], 
#         genes.metadata[names(sel.genes), , drop = FALSE]
#       )
#     )
#   }
# }

.filterGenesSparse <- function(
  counts,
  genes.metadata,
  gene.ID.column,
  cells.metadata,
  cell.type.column,
  filt.genes.cluster,
  min.mean.counts,
  filt.genes.cells,
  min.counts,
  min.cells,
  verbose
) {
  # duplicated genes in count matrix (and genes.metadata)
  dup.genes <- duplicated(rownames(counts))
  if (any(dup.genes)) {
    if (verbose) {
      message("      - Aggregating ", sum(dup.genes), " duplicated genes\n") 
    }
    ## this part will be changed
    counts <- rowsum(as.matrix(counts), rownames(counts))
  }
  genes.metadata <- genes.metadata[match(rownames(counts), 
                                         genes.metadata[, gene.ID.column]), , 
                                   drop = FALSE]
  # removing genes with no expression
  row.zero <- Matrix::rowSums(counts) > 0
  if (!all(row.zero)) {
    if (verbose) {
      message(paste("      - Removing", sum(!row.zero),
                    "genes without expression in any cell")) 
    }
    counts <- counts[row.zero, ]
    genes.metadata <- genes.metadata[genes.metadata[, gene.ID.column] %in%
                                       rownames(counts), , drop = FALSE]
  }
  dim.bef <- dim(counts)
  # filtering genes by expression thresholds: if both criteria are used, 
  # unexpected results can happen
  if (filt.genes.cells) {
    list.data <- .filterGenesByCells(
      counts = counts, 
      genes.metadata = genes.metadata, 
      min.counts = min.counts, 
      min.cells = min.cells 
    )
    counts <- list.data[[1]]
    genes.metadata <- list.data[[2]]
  }
  # if (filt.genes.cluster) {
  #   list.data <- .filterGenesByCluster(
  #     counts = counts, 
  #     genes.metadata = genes.metadata, 
  #     cells.metadata = cells.metadata, 
  #     cell.type.column = cell.type.column, 
  #     min.mean.counts = min.mean.counts
  #   )
  #   counts <- list.data[[1]]
  #   genes.metadata <- list.data[[2]]
  # }
  if (dim(counts)[1] == 0) {
    stop(paste("Resulting count matrix after filtering does not have entries"))
  }
  if (verbose) {
    message("      - Filtering features:")
    message(paste("         - Selected features:",  dim(counts)[1]))
    message(paste("         - Discarded features:", dim.bef[1] - dim(counts)[1]))  
  }
  genes.metadata <- genes.metadata[genes.metadata[, gene.ID.column] %in%
                                     rownames(counts), , drop = FALSE]
  
  return(list(counts, genes.metadata))
}

.filterGenesHDF5 <- function(
  counts,
  genes.metadata,
  gene.ID.column,
  cells.metadata,
  cell.type.column,
  filt.genes.cells,
  min.counts,
  min.cells,
  filt.genes.cluster,
  min.mean.counts,
  verbose
) { 
  if (verbose) {
    message("\n=== Processing data in HDF5 by blocks\n")
  }
  dup.genes <- duplicated(rownames(counts))
  if (any(dup.genes)) {
    if (verbose) {
      message("    - Aggregating ", sum(dup.genes), " duplicated genes\n") 
    }
    counts <- DelayedArray::rowsum(x = counts, group = factor(rownames(counts)))
    genes.metadata <- genes.metadata[match(
      x = rownames(counts), table = genes.metadata[, gene.ID.column]
    ), ]
  }
  # removing genes without any expression
  row.zero <- DelayedArray::rowSums(counts) > 0
  if (!all(row.zero)) {
    if (verbose) {
      message(paste("    - Removing", sum(!row.zero),
                    "genes without expression in any cell\n"))  
    }
    counts <- counts[row.zero, ]
    if (is.null(rownames(counts))) {
      genes.metadata <- genes.metadata[row.zero, , drop = FALSE]
    } else {
      genes.metadata <- genes.metadata[genes.metadata[, gene.ID.column] %in%
                                         rownames(counts), , drop = FALSE]  
    }
  }
  # filtered genes by cells
  ori.features <- nrow(counts)
  if (filt.genes.cells) {
    if (min.counts < 0 || min.cells < 0) {
      stop("min.counts and min.cells must be greater than or equal to zero")
    } else if (min.counts != 0 && min.cells != 0) {
      remove.genes <- DelayedArray::rowSums(counts > min.counts) >= min.cells
      counts <- counts[remove.genes, ]
      if (dim(counts)[1] == 0) {
        stop(paste("Resulting count matrix after filtering using min.genes =",
                   min.counts, "and min.cells =", min.cells,
                   "does not have entries"))
      }
      if (is.null(rownames(counts))) {
        genes.metadata <- genes.metadata[remove.genes, , drop = FALSE]
      } else {
        genes.metadata <- genes.metadata[genes.metadata[, gene.ID.column] %in%
                                           rownames(counts), , drop = FALSE]  
      }  
    }
  }
  # filter genes by cluster
  # if (filt.genes.cluster) {
  #   sum.cluster <- DelayedArray::rowsum(
  #     x = DelayedArray::t(counts), group = factor(cells.metadata[[cell.type.column]])
  #   )
  #   n.cluster <- data.frame(table(cells.metadata[[cell.type.column]]))
  #   mean.cluster <- sum.cluster[n.cluster$Var1,] / n.cluster$Freq
  #   
  #   ## using cutoffs
  #   sel.genes <- unlist(
  #     apply(
  #       X = mean.cluster >= min.mean.counts, 
  #       MARGIN = 2, 
  #       FUN = function(x) if(any(x)) return(TRUE)
  #     )
  #   )  
  #   if (any(sel.genes)) {
  #     counts <- counts[names(sel.genes), ]
  #     genes.metadata <- genes.metadata[names(sel.genes), ]
  #   }
  # }
  final.features <- nrow(counts)
  if (verbose) {
    message("\n    - Filtering features:")
    message(paste("       - Selected features:",  final.features))
    message(paste("       - Discarded features:", ori.features - final.features)) 
  }
  return(list(counts, genes.metadata))
}

.randomStr <- function() {
  a <- do.call(paste0, replicate(5, sample(LETTERS, 1, TRUE), FALSE))
  return(paste0("/", a, sprintf("%04d", sample(9999, 1, TRUE)), 
                sample(LETTERS, 1, TRUE)))
}

.loadSCData <- function(
  single.cell, 
  cell.ID.column,
  cell.type.column,
  gene.ID.column,
  name.dataset.h5,
  filt.genes.cluster = TRUE,
  min.mean.counts = 0,
  filt.genes.cells = TRUE,
  min.cells = 0, 
  min.counts = 0,
  file.backend = NULL,
  name.dataset.backend = NULL,
  compression.level = NULL,
  chunk.dims = NULL,
  block.processing = FALSE,
  verbose = TRUE
) {
  if (is.null(single.cell)) {
    stop(paste("Please, provide a 'single.cell' argument"))
  } else if (missing(cell.ID.column) || missing(gene.ID.column) || 
             is.null(cell.ID.column) || is.null(gene.ID.column)) {
    stop("'cell.ID.column' and 'gene.ID.column' arguments are needed. Please, see ",
         "?loadSCProfiles")
  }
  if (!is.null(file.backend)) {
    hdf5Params <- .checkHDF5parameters(
      file.backend = file.backend, 
      name.dataset.backend = name.dataset.backend, 
      compression.level = compression.level
    )
    name.dataset.backend <- hdf5Params[[1]]
    compression.level <- hdf5Params[[2]]
  }
  if (verbose) {
    message("=== Processing single-cell data")
  }
  if (is(single.cell, "SingleCellExperiment")) {
    # extract data (no filtering)
    list.data <- .extractDataFromSE(
      SEobject = single.cell,
      cell.ID.column = cell.ID.column,
      gene.ID.column = gene.ID.column,
      min.counts = min.counts,
      min.cells = min.cells
    )
  } else if (length(single.cell) == 0) {
    stop(paste("'single.cell' argument is empty"))
  } else if (length(single.cell) == 3 && !missing(name.dataset.h5)) {
    # from file --> hdf5 (needs dataset name)
    list.data <- list(
      .readCountsFile(
        counts.file = single.cell[[1]], 
        name.h5 = name.dataset.h5, 
        file.backend = file.backend,
        block.processing = block.processing
      ),
      .readTabFiles(single.cell[[2]]),
      .readTabFiles(single.cell[[3]])
    )
  } else if (length(single.cell) == 3) {
    # from files --> tsv, tsv.gz, mtx
    list.data <- list(
      .readCountsFile(single.cell[[1]]),
      .readTabFiles(single.cell[[2]]),
      .readTabFiles(single.cell[[3]])
    )
  } else {
    stop("Incorrect number of data elements given. Please, look at ", 
         "allowed data for in ?loadSCProfiles")
  }
  # use HDF5 backend and block.processing from both SCE object and files
  if (block.processing && is.null(file.backend)) {
    stop("block.processing is only compatible with HDF5 files used as back-end") 
  } else if (block.processing && !is.null(file.backend)) {
    if (!class(list.data[[1]]) %in% c("HDF5Matrix", "HDF5Array", 
                                      "DelayedArray", "DelayedMatrix")) {
      if (verbose) {
        message("=== Provided data are not stored as HDF5 file and ", 
                "'block.processing' has been set to TRUE, so data will be ", 
                "written in HDF5 file for block processing")
      }
      list.data[[1]] <- .useH5backend(
        counts = list.data[[1]], 
        file.backend = HDF5Array::getHDF5DumpFile(for.use = TRUE),
        compression.level = compression.level,
        group = HDF5Array::getHDF5DumpName(for.use = TRUE),
        # verbose = verbose
      ) 
    }
  } else if (!block.processing) {
    ## just in case the element provided is not a sparse matrix
    if (is(list.data[[1]], "matrix") | is(list.data[[1]], "data.frame")) {
      list.data[[1]] <- Matrix::Matrix(as.matrix(list.data[[1]]), sparse = TRUE)  
    } else if (is(list.data[[1]], "dgTMatrix")) {
      list.data[[1]] <- as(list.data[[1]], "dgCMatrix") 
    }
  }
  list.data <- .processData(
    counts = list.data[[1]],
    cells.metadata = list.data[[2]],
    cell.ID.column = cell.ID.column,
    cell.type.column = cell.type.column,
    genes.metadata = list.data[[3]],
    gene.ID.column = gene.ID.column,
    filt.genes.cluster = filt.genes.cluster,
    min.mean.counts = min.mean.counts,
    filt.genes.cells = filt.genes.cells,
    min.counts = min.counts,
    min.cells = min.cells,
    block.processing = block.processing,
    verbose = verbose
  )
  return(
    .createSCEObject(
      counts = list.data[[1]],
      cells.metadata = list.data[[2]],
      genes.metadata = list.data[[3]],
      file.backend = file.backend,
      name.dataset.backend = name.dataset.backend,
      compression.level = compression.level,
      chunk.dims = chunk.dims,
      block.processing = block.processing,
      verbose = verbose
    )
  )
}

.extractDataFromSE <- function(
    SEobject,
    cell.ID.column,
    gene.ID.column,
    min.counts = 0,
    min.cells = 0,
    new.data = TRUE
) {
  # extract cells.metadata
  cells.metadata <- SummarizedExperiment::colData(SEobject)
  if (any(dim(cells.metadata) == 0)) {
    stop("No data provided in colData slot. Cells metadata are needed. ",
         "Please, see ?loadSCProfiles")
  }
  if (!missing(cell.ID.column) && new.data) {
    # check if given IDs exist in cells.metadata. In cells.metadata are not
    # necessary because the data are provided from an SCE object
    .checkColumn(
      metadata = cells.metadata,
      ID.column = cell.ID.column,
      type.metadata = "cells.metadata",
      arg = "cell.ID.column"
    )
  }
  # extract count matrix
  if (length(SummarizedExperiment::assays(SEobject)) == 0) {
    stop("No count data in SingleCellExperiment object")
  } else if (length(SummarizedExperiment::assays(SEobject)) > 1) {
    warning("There is more than one assay, only the first will be used. ", 
            "Remember it must be raw data and not log-transformed data")
  }
  counts <- SummarizedExperiment::assays(SEobject)[[1]]
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("Count matrix must have rownames corresponding to features and ",  
         "colnames corresponding to cells")
  }
  # extract genes.metadata
  genes.metadata <- SummarizedExperiment::rowData(SEobject)
  if (!missing(gene.ID.column) && new.data) {
    if (any(dim(genes.metadata) == 0)) {
      stop("No data provided in rowData slot. Genes metadata are needed. ",
           "Please, see ?loadSCProfiles")
      # if (class(gene.ID.column) == "numeric") gene.ID.column <- "gene_names"
      # genes.metadata <- S4Vectors::DataFrame(gene.ID.column = rownames(counts))
    }
    # check if given IDs exist in genes.metadata. In cells.metadata are not
    # necessary because the data are provided from a SCE object
    .checkColumn(
      metadata = genes.metadata,
      ID.column = gene.ID.column,
      type.metadata = "genes.metadata",
      arg = "gene.ID.column"
    )
  }
  return(list(counts, cells.metadata, genes.metadata))
}

.loadSTData <- function(
  list.st.objects,
  cell.ID.column, 
  gene.ID.column,
  min.cells = 0, 
  min.counts = 0,
  n.slides = 1,
  verbose = TRUE
) {
  if (is(list.st.objects, "SpatialExperiment")) {
    list.st.objects <- list(list.st.objects)
    n.slides <- 1
  } else if (!is(list.st.objects, "list")) {
    stop(
      "`st.data` argument must be a SpatialExperiment object ", 
      "or a list of SpatialExperiment objects"
    )
  }
  if (verbose) {
    message(
      paste0(
        "=== ", length(list.st.objects), " SpatialExperiment objects provided\n"
      )
    )
  }
  ## check if there are names in the provided list
  if (!is.null(names(list.st.objects))) {
    if (length(unique(names(list.st.objects))) != length(names(list.st.objects))) {
      warning(
        "`st.data` list contains repeated names. Making them unique", 
        call. = FALSE, immediate. = TRUE
      )
      names.st.objs <- make.unique(names(list.st.objects))
    } else {
      names.st.objs <- names(list.st.objects)
    }
  } else {
    names.st.objs <- NULL
  }
  list.st.objects <- lapply(
    X = list.st.objects,
    FUN = function(obj) {
      obj.mod <- .processSTData(
        st.object = obj, 
        cell.ID.column = cell.ID.column, 
        gene.ID.column = gene.ID.column,
        min.cells = min.cells, 
        min.counts = min.counts,
        verbose = TRUE
      )
      genes <- rowData(obj.mod)[[gene.ID.column]]
      return(list(obj.mod, genes))
    }
  )
  # TODO: include messages with number of genes lost, etc.
  ## calculate genes shared in at least X slides (assuming genes are unique in each vector)
  if (n.slides == 0 | n.slides > length(list.st.objects)) {
    warning(
      paste0("`st.n.slides` parameter must be a valid number between 1 and the total", 
      " number of SpatialExperiment objects provided in `st.data`. 
      Setting st.n.slides = ", length(list.st.objects))
    )
    n.slides <- length(list.st.objects)
  }
  ## if n.slides == 1 & length(list.st.objects) == 1: only the first element, 
  # no subset
  if (n.slides == 1 & length(list.st.objects) == 1) {
    return(
      list(list.st.objects[[1]][[1]]) %>% setNames(names.st.objs)
    )
  ## if n.slides == 1 but length(list.st.objects) != 1: no subset, but lapply
  # to make it a list of n elements
  } else if (n.slides == 1) {
    return(
      lapply(X = list.st.objects, FUN = function(x) x[[1]]) %>% 
        setNames(names.st.objs)
    )
  ## if n.slides != 1 & length(list.st.objects) != 1: subset, vectorized
  } else {
    all.genes <- sapply(X = list.st.objects, FUN = function(x) x[[2]])
    num.genes <- table(unlist(all.genes))
    inter.genes <- names(num.genes[num.genes >= n.slides])

    return(
      lapply(
        list.st.objects, function(x) {
          ## genes not present in inter.genes: they will be 0s
          genes.out <- inter.genes[!inter.genes %in% x[[2]]]
          ## now, inter.genes.subset will be used to filter the original matrix
          inter.genes.subset <- inter.genes[inter.genes %in% x[[2]]]
          x[[1]] <- x[[1]][inter.genes.subset, ]
          ## and now we can generate rows for those genes not detected in this slide
          if (identical(genes.out, character(0))) {
            return(x[[1]])
          } else {
            zero.genes <- matrix(
              0L, nrow = length(genes.out), ncol = ncol(x[[1]]), 
              dimnames = list(genes.out, NULL)
            )
            counts <- rbind(
              SummarizedExperiment::assays(x[[1]])[[1]], zero.genes
            )[inter.genes, ]
            genes.metadata <- S4Vectors::DataFrame(rownames(counts)) %>% 
              setNames(colnames(SummarizedExperiment::rowData(x[[1]])))
            
            return(
              SpatialExperiment::SpatialExperiment(
                assays = rbind(
                  SummarizedExperiment::assays(x[[1]])[[1]], zero.genes
                )[inter.genes, ],
                colData = x[[1]]@colData,
                rowData = S4Vectors::DataFrame(rownames(counts)) %>% 
                  setNames(colnames(SummarizedExperiment::rowData(x[[1]]))),
                spatialCoords = SpatialExperiment::spatialCoords(x[[1]])
              )
            )
          }
        }
      ) %>% setNames(names.st.objs)
    )
  }
}


## set of functions to load ST data
.processSTData <- function(
    st.object, 
    cell.ID.column, 
    gene.ID.column,
    min.cells = 0, 
    min.counts = 0,
    verbose = TRUE
) {
  if (is.null(st.object)) {
    stop(paste("Please, provide a 'st.object' argument"))
  } else if (missing(cell.ID.column) || missing(gene.ID.column) || 
             is.null(cell.ID.column) || is.null(gene.ID.column)) {
    stop("'cell.ID.column' and 'gene.ID.column' arguments are needed. Please, see ",
         "?createSpatialDDLSobject or ?loadSTProfiles")
  }
  if (!is(st.object, "SpatialExperiment")) {
    stop(
      "Only SpatialExperiment objects are supported. See",
      " ?createSpatialDDLSobject or ?loadSTProfiles for more details"
    ) 
  } 
  # extract data
  # this function works well with SpatialExperiment too
  list.data <- .extractDataFromSE(
    SEobject = st.object,
    cell.ID.column = cell.ID.column,
    gene.ID.column = gene.ID.column,
    min.counts = min.counts,
    min.cells = min.cells
  )
  ## just in case the element provided is not a sparse matrix
  if (is(list.data[[1]], "matrix") | is(list.data[[1]], "data.frame")) {
    list.data[[1]] <- Matrix::Matrix(as.matrix(list.data[[1]]), sparse = TRUE)  
  } else if (is(list.data[[1]], "dgTMatrix")) {
    list.data[[1]] <- as(list.data[[1]], "dgCMatrix") 
  }
  if (verbose) {
    message("   === Processing spatial transcriptomics data")
  }
  ## filtering genes: setting not used arguments manually
  list.data <- .processData(
    counts = list.data[[1]],
    cells.metadata = list.data[[2]],
    cell.ID.column = cell.ID.column,
    cell.type.column = NULL,
    genes.metadata = list.data[[3]],
    gene.ID.column = gene.ID.column,
    filt.genes.cluster = FALSE,
    min.mean.counts = NULL,
    filt.genes.cells = TRUE,
    min.counts = min.counts,
    min.cells = min.cells,
    block.processing = FALSE,
    file.backend = NULL,
    verbose = verbose
  )
  # ## modify original st.object
  # st.object@assays <- SummarizedExperiment::Assays(list.data[[1]])
  # st.object@colData <- list.data[[2]]
  # # TODO: error
  # rowData(st.object) <- list.data[[3]]
  # message(class(list.data[[1]]))
  
  ## just for a better visualization
  if (verbose) message("")
  
  return(
    SpatialExperiment::SpatialExperiment(
      assays = list.data[[1]],
      colData = list.data[[2]],
      rowData = list.data[[3]],
      spatialCoords = SpatialExperiment::spatialCoords(st.object)
    )
  )
}

.filterGenesByCluster <- function(
  sce.obj,
  cell.type.column, 
  min.mean.counts,
  n.genes.per.cluster,
  top.n.genes,
  verbose
) {
  
  if (is.null(names(sce.obj@assays@data))) 
    names(sce.obj@assays@data) <- "counts"
  
  mean.cluster.counts <- .aggregate.Matrix.sparse(
    x = Matrix::t(sce.obj@assays@data$counts), 
    groupings = list(sce.obj[[cell.type.column]]), 
    fun = "mean"
  )
  sel.genes <- unlist(
    apply(
      X = mean.cluster.counts > min.mean.counts, 
      MARGIN = 2, 
      FUN = function(x) if(any(x)) return(TRUE)
    )
  )
  sce.obj.norm <- scuttle::computeLibraryFactors(sce.obj[names(sel.genes), ])
  sce.obj.norm <- scuttle::logNormCounts(sce.obj.norm)
  ## mean values
  mean.cluster.log <- .aggregate.Matrix.sparse(
    x = Matrix::t(sce.obj.norm@assays@data$logcounts), 
    groupings = list(sce.obj.norm@colData[[cell.type.column]]), 
    fun = "mean"
  )
  means.all <- colMeans(as.matrix(mean.cluster.log))
  logFCs.cluster <- apply(mean.cluster.log, 1, \(x) x - means.all)
  list.cluster.FC <- lapply(
    as.list(as.data.frame(logFCs.cluster)), 
    \(x) x %>% setNames(rownames(logFCs.cluster))
  )
  ranked.logFC.top <- lapply(
    list.cluster.FC, 
    \(x) {
      x.f <- x[x > 0.5] ## only if logFC > 0.5
      x.f[order(x.f, decreasing = T)] %>% head(n.genes.per.cluster) %>% names()
    }
  ) %>% unlist() %>% unique()
  final.genes <- intersect(ranked.logFC.top, names(sel.genes))
  
  if (verbose) 
    message(
      "\n=== Number of genes after filtering based on logFC: ", length(final.genes)
    )
  
  if (length(final.genes) > top.n.genes) {
    if (verbose) 
      message(
        "\n=== As the number of resulting genes is greater than ",
        "the top.n.genes parameter. Using only ", 
        top.n.genes, " according to gene variance"
      )
    dec.sc.obj.ln <- scran::modelGeneVar(sce.obj.norm[final.genes, ])
    final.genes <- as.data.frame(dec.sc.obj.ln) %>% 
      arrange(desc(abs(.data[["bio"]]))) %>% head(top.n.genes) %>% rownames()
  }
  
  return(final.genes)
}

################################################################################
######################### Create a SpatialDDLS object ##########################
################################################################################

#' Create a \code{\linkS4class{SpatialDDLS}} object
#'
#' Create a \code{\linkS4class{SpatialDDLS}} object by providing single-cell
#' RNA-seq data. Additionally, spatial transcriptomics data contained in
#' \code{\linkS4class{SpatialDDLS}} objects can also be provided. It is
#' recommended to provide both types of data to only use genes shared between 
#' both experiments. 
#' 
#' \strong{Filtering genes}
#' 
#' In order to reduce the number of dimensions used for subsequent steps, 
#' \code{createSpatialDDLSobject} implements different strategies aimed at 
#' removing useless genes for deconvolution: \itemize{ \item Filtering at the 
#' cell level: genes less expressed than a determined cutoff in N cells are
#' removed. See \code{sc.min.cells}/\code{st.min.cells} and 
#' \code{sc.min.counts}/\code{st.min.cells} parameters. \item Filtering at the 
#' cluster level (only for scRNA-seq data): if 
#' \code{sc.filt.genes.cluster == TRUE}, \code{createSpatialDDLSobject} sets a 
#' cutoff of non-zero average counts per 
#' cluster (\code{sc.min.mean.counts} parameter) and take only the 
#' \code{sc.n.genes.per.cluster} genes with the highest logFC per cluster. 
#' LogFCs are calculated using normalized logCPM of each cluster with respect to 
#' the average in the whole dataset). Finally, if 
#' the number of remaining genes is greater than \code{top.n.genes}, genes are 
#' ranked based on variance and the \code{top.n.genes} most variable genes are 
#' used for downstream analyses.}
#'
#' \strong{Single-cell RNA-seq data}
#'
#' Single-cell RNA-seq data can be provided from files (formats allowed: tsv,
#' tsv.gz, mtx (sparse matrix) and hdf5) or a
#' \code{\linkS4class{SingleCellExperiment}} object. Data will be stored in the
#' \code{single.cell.real} slot, and must consist of three pieces of
#' information: \itemize{ \item Single-cell counts: genes as rows and cells as
#' columns. \item Cells metadata: annotations (columns) for each cell (rows).
#' \item Genes metadata: annotations (columns) for each gene (rows). } If data
#' are provided from files, \code{single.cell.real} argument must be a vector of
#' three elements ordered so that the first file corresponds to the count
#' matrix, the second to the cells metadata, and the last to the genes metadata.
#' On the other hand, if data are provided as a
#' \code{\linkS4class{SingleCellExperiment}} object, it must contain single-cell
#' counts in \code{assay}, cells metadata in \code{colData}, and genes metadata
#' in \code{rowData}. Data must be provided without any transformation (e.g.
#' log-transformation), raw counts are preferred.
#'
#' \strong{Spatial transcriptomics data}
#'
#' It must be a \code{\linkS4class{SpatialExperiment}} object (or a list of them
#' if more than one slide is going to be deconvoluted) containing the same
#' information as the single-cell RNA-seq data: the count matrix, spots
#' metadata, and genes metadata. Please, make sure the gene identifiers used the
#' spatial and single-cell transcriptomics data are consistent.
#'
#' @param sc.data Single-cell RNA-seq profiles to be used as reference. If data
#'   are provided from files, \code{single.cell.real} must be a vector of three
#'   elements: single-cell counts, cells metadata and genes metadata. On the
#'   other hand, If data are provided from a
#'   \code{\linkS4class{SingleCellExperiment}} object, single-cell counts must
#'   be present in the \code{assay} slot, cells metadata in the \code{colData}
#'   slot, and genes metadata in the \code{rowData} slot.
#' @param sc.cell.ID.column Name or number of the column in cells metadata
#'   corresponding to cell names in expression matrix (single-cell RNA-seq
#'   data).
#' @param sc.cell.type.column Name or column number corresponding to cell types
#'   in cells metadata.
#' @param sc.gene.ID.column Name or number of the column in genes metadata
#'   corresponding to the names used for features/genes (single-cell RNA-seq
#'   data).
#' @param st.data Spatial transcriptomics datasets to be deconvoluted. It can be
#'   a single \code{\linkS4class{SpatialExperiment}} object or a list of them.
#' @param st.spot.ID.column Name or number of the column in spots metadata
#'   corresponding to spot names in expression matrix (spatial transcriptomics
#'   data).
#' @param st.gene.ID.column Name or number of the column in the genes metadata
#'   corresponding to the names used for features/genes (spatial transcriptomics
#'   data).
#' @param filter.mt.genes Regular expression matching mitochondrial genes to 
#'   be ruled out (\code{^mt-} by default). If \code{NULL}, no filtering is 
#'   performed. 
#' @param sc.filt.genes.cluster Whether to filter single-cell RNA-seq genes
#'   according to a minimum threshold of non-zero average counts per cell type
#'   (\code{sc.min.mean.counts}). \code{TRUE} by default. 
#' @param sc.min.mean.counts Minimum non-zero average counts per cluster to
#'   filter genes. 1 by default. 
#' @param sc.n.genes.per.cluster Top n genes with the highest logFC per cluster
#'   (300 by default). See Details section for more details. 
#' @param top.n.genes Maximum number of genes used for downstream steps (2000 
#'   by default). In case the number of genes after filtering is greater than 
#'   \code{top.n.genes}, these genes will be set according to 
#'   variability across the whole single-cell dataset. 
#' @param sc.min.counts Minimum gene counts to filter (1 by default; single-cell
#'   RNA-seq data).
#' @param sc.min.cells Minimum of cells with more than \code{min.counts} (1 by
#'   default; single-cell RNA-seq data).
#' @param st.min.counts Minimum gene counts to filter (1 by default; spatial
#'   transcriptomics data).
#' @param st.min.spots Minimum of cells with more than \code{min.counts} (1 by
#'   default; spatial transcriptomics data).
#' @param st.n.slides Minimum number of slides
#'   (\code{\linkS4class{SpatialExperiment}} objects) in which a gene has to be
#'   expressed in order to keep it. This parameter is applicable only when
#'   multiple \code{\linkS4class{SpatialExperiment}} objects are provided. Genes
#'   not present in at least \code{st.n.slides} will be discarded. If no
#'   filtering is desired, set \code{st.n.slides = 1}.
#' @param shared.genes If set to \code{TRUE}, only genes present in both the
#'   single-cell and spatial transcriptomics data will be retained for further
#'   processing (\code{TRUE} by default).
#' @param sc.name.dataset.h5 Name of the data set if HDF5 file is provided for
#'   single-cell RNA-seq data.
#' @param sc.file.backend Valid file path where to store the loaded for
#'   single-cell RNA-seq data as HDF5 file. If provided, data are stored in a
#'   HDF5 file as back-end using the \pkg{DelayedArray} and \pkg{HDF5Array}
#'   packages instead of being loaded into RAM. This is suitable for situations
#'   where you have large amounts of data that cannot be stored in memory. Note
#'   that operations on these data will be performed by blocks (i.e subsets of
#'   determined size), which may result in longer execution times. \code{NULL}
#'   by default.
#' @param sc.name.dataset.backend Name of the HDF5 file dataset to be used. Note
#'   that it cannot exist. If \code{NULL} (by default), a random dataset name
#'   will be generated.
#' @param sc.compression.level The compression level used if
#'   \code{sc.file.backend} is provided. It is an integer value between 0 (no
#'   compression) and 9 (highest and slowest compression). See
#'   \code{?\link[HDF5Array]{getHDF5DumpCompressionLevel}} from the
#'   \pkg{HDF5Array} package for more information.
#' @param sc.chunk.dims Specifies dimensions that HDF5 chunk will have. If
#'   \code{NULL}, the default value is a vector of two items: the number of
#'   genes considered by \code{\linkS4class{SpatialDDLS}} object during the
#'   simulation, and only one sample in order to increase read times in the
#'   following steps. A larger number of columns written in each chunk may lead
#'   to longer read times.
#' @param sc.block.processing Boolean indicating whether single-cell RNA-seq
#'   data should be treated as blocks (only if data are provided as HDF5 file).
#'   \code{FALSE} by default. Note that using this functionality is suitable for
#'   cases where it is not possible to load data into RAM and therefore
#'   execution times will be longer.
#' @param verbose Show informative messages during the execution (\code{TRUE} by
#'   default).
#' @param project Name of the project for \code{\linkS4class{SpatialDDLS}}
#'   object.
#'
#' @return A \code{\linkS4class{SpatialDDLS}} object with the single-cell
#'   RNA-seq data provided loaded into the \code{single.cell.real} slot as a
#'   \code{\linkS4class{SingleCellExperiment}} object. If spatial
#'   transcriptomics data are provided, they will be loaded into the
#'   \code{spatial.experiments} slot.
#'
#' @export
#'
#' @seealso \code{\link{estimateZinbwaveParams}} \code{\link{genMixedCellProp}}
#'
#' @examples
#' set.seed(123)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(
#'     counts = matrix(
#'       rpois(100, lambda = 5), nrow = 40, ncol = 30,
#'       dimnames = list(paste0("Gene", seq(40)), paste0("RHC", seq(30)))
#'     )
#'   ),
#'   colData = data.frame(
#'     Cell_ID = paste0("RHC", seq(30)),
#'     Cell_Type = sample(x = paste0("CellType", seq(4)), size = 30,
#'                        replace = TRUE)
#'   ),
#'   rowData = data.frame(
#'     Gene_ID = paste0("Gene", seq(40))
#'   )
#' )
#' counts <- matrix(
#'   rpois(30, lambda = 5), ncol = 6,
#'   dimnames = list(paste0("Gene", 1:5), paste0("Spot", 1:6))
#' )
#' coordinates <- matrix(
#'   c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3), ncol = 2
#' )
#' ste <- SpatialExperiment::SpatialExperiment(
#'   assays = list(counts = as.matrix(counts)),
#'   rowData = data.frame(Gene_ID = paste0("Gene", 1:5)),
#'   colData = data.frame(Cell_ID = paste0("Spot", 1:6)),
#'   spatialCoords = coordinates
#' )
#'
#' SDDLS <- createSpatialDDLSobject(
#'   sc.data = sce,
#'   sc.cell.ID.column = "Cell_ID",
#'   sc.gene.ID.column = "Gene_ID",
#'   st.data = ste,
#'   st.spot.ID.column = "Cell_ID",
#'   st.gene.ID.column = "Gene_ID",
#'   project = "Simul_example",
#'   sc.filt.genes.cluster = FALSE
#' )
#'   
createSpatialDDLSobject <- function(
    sc.data,
    sc.cell.ID.column,
    sc.cell.type.column,
    sc.gene.ID.column,
    st.data,
    st.spot.ID.column,
    st.gene.ID.column,
    filter.mt.genes = "^mt-",
    sc.filt.genes.cluster = TRUE,
    sc.min.mean.counts = 1, 
    sc.n.genes.per.cluster = 300,
    top.n.genes = 2000,
    sc.min.counts = 1,
    sc.min.cells = 1,
    st.min.counts = 1,
    st.min.spots = 1,
    st.n.slides = 3,
    shared.genes = TRUE,
    sc.name.dataset.h5 = NULL, 
    sc.file.backend = NULL,
    sc.name.dataset.backend = NULL,
    sc.compression.level = NULL,
    sc.chunk.dims = NULL,
    sc.block.processing = FALSE,
    verbose = TRUE,
    project = "SpatialDDLS-Proj"
) {
  if (missing(sc.cell.type.column)) sc.cell.type.column <- NULL
  # in case filtering according to expression in each cluster is used
  if (sc.filt.genes.cluster & (is.null(sc.cell.type.column) | missing(sc.cell.type.column))) {
    stop("sc.cell.type.column must be provided")
  } 
  ## spatial transcriptomics profiles
  if (!missing(st.data)) {
    spatial.experiments <- .loadSTData(
      list.st.objects = st.data,
      cell.ID.column = st.spot.ID.column,
      gene.ID.column = st.gene.ID.column,
      min.cells = st.min.spots,
      min.counts = st.min.counts,
      n.slides = st.n.slides,
      verbose = verbose
    ) # %>% lapply(X = ., function(x) x[[1]])
  } else {
    spatial.experiments <- NULL
    if (verbose) message("=== Spatial transcriptomics data not provided")
  }
  ## single-cell profiles
  single.cell.real <- .loadSCData(
    single.cell = sc.data,
    cell.ID.column = sc.cell.ID.column,
    gene.ID.column = sc.gene.ID.column,
    cell.type.column = sc.cell.type.column,
    name.dataset.h5 = sc.name.dataset.h5,
    filt.genes.cluster = FALSE,
    min.mean.counts = sc.min.mean.counts,
    filt.genes.cells = TRUE,
    min.cells = sc.min.cells,
    min.counts = sc.min.counts,
    file.backend = sc.file.backend,
    name.dataset.backend = sc.name.dataset.backend,
    compression.level = sc.compression.level,
    chunk.dims = sc.chunk.dims,
    block.processing = sc.block.processing,
    verbose = verbose
  )
## intersection between datasets
  if (!missing(st.data)) {
    if (shared.genes) {
      ## slide with more genes
      pos.max <- which.max(sapply(spatial.experiments, nrow))
      inter.genes <- intersect(
        rownames(single.cell.real), rownames(spatial.experiments[[pos.max]])
      )
      if (verbose) {
        message(
          "\n=== Number of shared genes between single-cell and spatial transcriptomics datasets: ", 
          length(inter.genes) 
        )
        message(
          "    - Original # genes in single-cell data: ", nrow(single.cell.real) 
        )
        message(
          "    - Original # genes in ST data (object with the greatest # genes): ", 
          nrow(spatial.experiments[[pos.max]]) 
        )
      }
      single.cell.real <- single.cell.real[inter.genes, ]
      spatial.experiments <- lapply(
        X = spatial.experiments, 
        FUN = function(obj) obj[inter.genes[inter.genes %in% rownames(obj)], ]
      )
    }
  }
  ## take out mitochondrial genes
  if (!is.null(filter.mt.genes)) {
    mt.genes.sc <- grepl(
      pattern = filter.mt.genes, x = rownames(single.cell.real), 
      ignore.case = TRUE
    )
    if (sum(mt.genes.sc) > 0) {
      if (verbose) 
        message("\n=== Number of removed mitochondrial genes: ", sum(mt.genes.sc))
      
      single.cell.real <- single.cell.real[!mt.genes.sc, ]
      
      mt.genes.st <- grepl(
        pattern = filter.mt.genes, x = rownames(spatial.experiments[[1]]), 
        ignore.case = TRUE
      )
      spatial.experiments <- lapply(spatial.experiments, \(st) st[!mt.genes.st, ])
    } else {
      if (verbose) 
        message(
          "\n=== No mitochondrial genes were found by using ", 
          filter.mt.genes, " as regrex"
        ) 
    }
  }
  if (sc.filt.genes.cluster) {
    final.genes <- .filterGenesByCluster(
      sce.obj = single.cell.real,
      cell.type.column = sc.cell.type.column, 
      min.mean.counts = sc.min.mean.counts,
      n.genes.per.cluster = sc.n.genes.per.cluster,
      top.n.genes = top.n.genes,
      verbose = verbose
    )  
    spatial.experiments <- lapply(spatial.experiments, \(st) st[final.genes, ])
    single.cell.real <- single.cell.real[final.genes, ]
  }
  
  ## messages
  if (verbose) 
    message(
      "\n=== Final number of dimensions for further analyses: ", 
      nrow(single.cell.real)
    )
  
  return(
    new(
      Class = "SpatialDDLS",
      single.cell.real = single.cell.real,
      spatial.experiments = spatial.experiments,
      project = project,
      version = packageVersion(pkg = "SpatialDDLS")
    )
  )
}


################################################################################
############################### Load spatial data ##############################
################################################################################

#' Loads spatial transcriptomics data into a SpatialDDLS object
#'
#' This function loads a \code{\linkS4class{SpatialExperiment}} object (or a
#' list with several \code{\linkS4class{SpatialExperiment}} objects) into a
#' \code{\linkS4class{SpatialDDLS}} object.
#'
#' It is recommended to perform this step when creating the
#' \code{\linkS4class{SpatialDDLS}} object using the
#' \code{\link{createSpatialDDLSobject}} function in order to only keep genes
#' shared between the spatial transcriptomics and the single-cell
#' transcriptomics data used as reference. In addition, please, make sure the
#' gene identifiers used the spatial and single-cell transcriptomics data are
#' consistent.
#'
#' @param object A \code{\linkS4class{SpatialDDLS}} object.
#' @param st.data A \code{\linkS4class{SpatialExperiment}} object (or a list
#'   with several \code{\linkS4class{SpatialExperiment}} objects) to be
#'   deconvoluted.
#' @param st.spot.ID.column Name or number of the column in spots metadata
#'   corresponding to spot names in the expression matrix.
#' @param st.gene.ID.column Name or number of the column in genes metadata
#'   corresponding to names used for features/genes.
#' @param st.min.counts Minimum gene counts to filter (0 by default).
#' @param st.min.spots Minimum of spots with more than \code{min.counts} (0 by
#'   default).
#' @param st.n.slides Minimum number of slides
#'   (\code{\linkS4class{SpatialExperiment}} objects) in which a gene has to be
#'   expressed in order to keep it. This parameter is applicable only when
#'   multiple \code{\linkS4class{SpatialExperiment}} objects are provided. Genes
#'   not present in at least \code{st.n.slides} will be discarded. If no
#'   filtering is desired, set \code{st.n.slides = 1}.
#' @param verbose Show informative messages during execution (\code{TRUE} by
#'   default).
#'
#' @return A \code{\linkS4class{SpatialDDLS}} object with the provided spatial
#'   trainscriptomics data loaded into the \code{spatial.experiments} slot.
#'
#' @export
#'
#' @seealso \code{\link{createSpatialDDLSobject}} \code{\link{trainDeconvModel}}
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(
#'     counts = matrix(
#'       rpois(100, lambda = 5), nrow = 40, ncol = 30,
#'       dimnames = list(paste0("Gene", seq(40)), paste0("RHC", seq(30)))
#'     )
#'   ),
#'   colData = data.frame(
#'     Cell_ID = paste0("RHC", seq(30)),
#'     Cell_Type = sample(x = paste0("CellType", seq(4)), size = 30,
#'                        replace = TRUE)
#'   ),
#'   rowData = data.frame(
#'     Gene_ID = paste0("Gene", seq(40))
#'   )
#' )
#' SDDLS <- createSpatialDDLSobject(
#'   sc.data = sce,
#'   sc.cell.ID.column = "Cell_ID",
#'   sc.gene.ID.column = "Gene_ID",
#'   sc.filt.genes.cluster = FALSE
#' )
#' 
#' ## simulating a SpatialExperiment object
#' counts <- matrix(rpois(30, lambda = 5), ncol = 6)
#' rownames(counts) <- paste0("Gene", 1:5)
#' colnames(counts) <- paste0("Spot", 1:6)
#' coordinates <- matrix(
#'   c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3), ncol = 2
#' )
#' ste <- SpatialExperiment::SpatialExperiment(
#'   assays = list(counts = as.matrix(counts)),
#'   rowData = data.frame(Gene_ID = paste0("Gene", 1:5)),
#'   colData = data.frame(Cell_ID = paste0("Spot", 1:6)),
#'   spatialCoords = coordinates
#' )
#'
#' ## previous SpatialDDLS object
#' SDDLS <- loadSTProfiles(
#'   object = SDDLS,
#'   st.data = ste,
#'   st.spot.ID.column = "Cell_ID",
#'   st.gene.ID.column = "Gene_ID"
#' )
#' }
#'   
loadSTProfiles <- function(
    object,
    st.data,
    st.spot.ID.column,
    st.gene.ID.column,
    st.min.counts = 0,
    st.min.spots = 0,
    st.n.slides = 3,
    verbose = TRUE
) {
  if (!is(object, "SpatialDDLS")) {
    stop("The object provided is not of SpatialDDLS class")
  } 
  spatial.experiments <- .loadSTData(
    list.st.objects = st.data,
    cell.ID.column = st.spot.ID.column,
    gene.ID.column = st.gene.ID.column,
    min.cells = st.min.spots,
    min.counts = st.min.counts,
    n.slides = st.n.slides,
    verbose = verbose
  )
  spatial.experiments(object) <- spatial.experiments
  
  return(object)
}
