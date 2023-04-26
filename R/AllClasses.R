#' @importFrom methods setClass setOldClass setClassUnion
#' @importFrom utils packageVersion
#' @import SingleCellExperiment SummarizedExperiment keras
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom stats predict
NULL

setOldClass(Classes = 'package_version')
# TODO: error with setOldClass, check what is going on
setClass("keras_training_history") 
setOldClass(Classes = "keras.engine.sequential.Sequential")

setClassUnion(name = "MatrixOrNULL", members = c("matrix", "NULL"))
setClassUnion(name = "ListOrNULL", members = c("list", "NULL"))
setClassUnion(name = "ListNumericOrNULL", members = c("list", "numeric", "NULL"))
setClassUnion(name = "CharacterOrNULL", members = c("character", "NULL"))
setClassUnion(name = "SingleCellExperimentOrNULL", 
              members = c("SingleCellExperiment", "NULL"))
setClassUnion(name = "KerasOrList", 
              members = c("keras.engine.sequential.Sequential", "list"))
setClassUnion(name = "KerasTrainOrNULL", 
              members = c("keras_training_history", "NULL"))

################################################################################
######################## Wrapper class for ZinbModel ###########################
################################################################################

#' The Class ZinbParametersModel
#'
#' The ZinbParametersModel class is a wrapper class for the
#' \code{\linkS4class{ZinbModel}} class from the \pkg{zinbwave} package.
#'
#' This wrapper class contains the \code{zinbwave.model} slot, which holds a
#' valid \code{\linkS4class{ZinbModel}} object.
#'
#' @slot zinbwave.model A valid \code{\linkS4class{ZinbModel}} object.
#'
#' @references Risso, D., Perraudeau, F., Gribkova, S. et al. (2018). A general
#'   and flexible method for signal extraction from single-cell RNA-seq data.
#'   Nat Commun 9, 284. doi: \doi{10.1038/s41467-017-02554-5}.
#'
#' @export ZinbParametersModel
#'   
ZinbParametersModel <- setClass(
  Class = "ZinbParametersModel",
  slots = c(zinbwave.model = "ANY")
)

setMethod(
  f = "initialize", 
  signature = "ZinbParametersModel",
  definition = function(
    .Object,
    zinbwave.model = NULL
  ) {
    .Object@zinbwave.model <- zinbwave.model
    return(.Object)
  }
)

setValidity(
  Class = "ZinbParametersModel",
  method = function(object) {
    if (object@zinbwave.model != "ZinbModel") {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
)

setMethod(
  f = "show",
  signature = "ZinbParametersModel",
  definition = function(object) {
    if (is.null(object@zinbwave.model)) {
      cat("ZinbParametersModel object empty")
    } else {
      # TODO: check if it works out
      .zinbModelShow(object@zinbwave.model) 
    }
  }
)

setClassUnion(
  name = "ZinbParametersModelOrNULL", members = c("ZinbParametersModel", "NULL")
)

################################################################################
############################# PropCellTypes class ##############################
################################################################################

#' The PropCellTypes Class
#'
#' The \code{\linkS4class{PropCellTypes}} class is a data storage class which
#' contains information related to the cell type composition matrices used to
#' simulate mixed transcriptional profiles. This matrix is stored in the
#' \code{prob.matrix} slot while the rest contains additional information
#' generated during the process and required in subsequent steps.
#'
#' See \code{?\link{genMixedCellProp}} function for information about how cell
#' type composition matrices are generated. Plots about the cell type proportion
#' distributions can be accessed using the \code{\link{showProbPlot}} function
#' (see \code{?\link{showProbPlot}} for more details).
#'
#' @slot prob.matrix Matrix of cell type proportions to simulate mixed
#'   transcriptional profiles.
#' @slot cell.names Matrix containing the name of cells to be used to generate
#'   the simulated mixed transcriptional profiles.
#' @slot set.list List of cells ordered according to the cell type they belong
#'   to.
#' @slot set Vector containing cell names present in the object.
#' @slot method Vector containing the method by which cell type proportions were
#'   generated.
#' @slot plots Plots showing the distribution of the cell type proportions. See
#'   \code{?\link{showProbPlot}} for more details.
#' @slot type.data Character indicating the type of data contained: train or
#'   test.
#'
#' @export PropCellTypes
#'   
PropCellTypes <- setClass(
  Class = "PropCellTypes",
  slots = c(
    prob.matrix = "MatrixOrNULL",
    cell.names = "MatrixOrNULL",
    set.list = "ListOrNULL",
    set = "CharacterOrNULL",
    method = "CharacterOrNULL",
    plots = "ListOrNULL",
    type.data = "CharacterOrNULL"
  )
)

setMethod(
  f = "initialize", 
  signature = "PropCellTypes",
  definition = function(
    .Object,
    prob.matrix = NULL,
    cell.names = NULL,
    set.list = NULL,
    set = NULL,
    method = NULL,
    plots = NULL,
    type.data = NULL
  ) {
    .Object@prob.matrix <- prob.matrix
    .Object@cell.names <- cell.names
    .Object@set.list <- set.list
    .Object@set <- set
    .Object@method <- method
    .Object@plots <- plots
    .Object@type.data <- type.data
    return(.Object)
  }
)

setValidity(
  Class = "PropCellTypes",
  method = function(object) {
    if (all(object@type.data != c("train", "test"))) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
)

setMethod(
  f = "show",
  signature = "PropCellTypes",
  definition = function(object) {
    # cat("An object of class", class(object), "\n")
    if (is.null(object@prob.matrix)) {
      cat("PropCellTypes object empty")
    } else {
      cat(paste0("  Cell type matrix for ", object@type.data, "data: "))
      cat(paste(dim(object@prob.matrix), 
                c("bulk samples and", "cell types"), collapse = " "))
    }
  }
)

################################################################################
############################ DeconvDLModel class ###############################
################################################################################

#' The DeconvDLModel Class
#'
#' The \code{\linkS4class{DeconvDLModel}} object stores all the information
#' related to deep neural network models. It consists of the trained model, the
#' training history, and the predictions on test data. After running
#' \code{\link{calculateEvalMetrics}}, it is possible to find the performance
#' evaluation of the model on test data (see \code{?\link{calculateEvalMetrics}}
#' for details).
#'
#' The steps related to Deep Learning are carried out using the \pkg{keras} and
#' \pkg{tensorflow} packages, which use the R6 classes system. If you want to
#' save the \code{\linkS4class{DeconvDLModel}} object as an RDS file,
#' \pkg{SpatialDDLS} provides a \code{saveRDS} generic function that transforms
#' the R6 object containing the trained model into a native valid R object.
#' Specifically, the model is converted into a list with the architecture of the
#' network and the weights learned during training. That is the minimum
#' information needed to use the model as predictor. If you want to keep the
#' optimizer state, see \code{?\link{saveTrainedModelAsH5}}. If you want to
#' store either the \code{\linkS4class{DeconvDLModel}} or the
#' \code{\linkS4class{SpatialDDLS}} objects on disk as RDA files, see
#' \code{?\link{preparingToSave}}.
#'
#' @slot model Trained deep neural network. This slot can contain an R6
#'   \code{keras.engine.sequential.Sequential} object or a list with two
#'   elements: the architecture of the model and the resulting weights after
#'   training.
#' @slot training.history List with the evolution of the selected metrics during
#'   training.
#' @slot test.metrics Performance of the model on test data.
#' @slot test.pred Predicted cell type proportions on test data.
#' @slot cell.types Vector with cell types considered by the model.
#' @slot features Vector with features (genes) considered by the model. These
#'   features will be used for subsequent predictions.
#' @slot test.deconv.metrics Performance of the model on test data by cell type.
#'   This slot is generated by the \code{\link{calculateEvalMetrics}} function
#'   (see \code{?\link{calculateEvalMetrics}} for more details).
#'
#' @export DeconvDLModel
#'   
DeconvDLModel <- setClass(
  Class = "DeconvDLModel",
  slots = c(
    model = "KerasOrList",
    training.history = "KerasTrainOrNULL",
    test.metrics = "ListNumericOrNULL",
    test.pred = "MatrixOrNULL",
    cell.types = "character",
    features = "character",
    test.deconv.metrics = "ListOrNULL"
  )
)

setMethod(
  f = "initialize", 
  signature = "DeconvDLModel",
  definition = function(
    .Object,
    model = list(),
    training.history = NULL,
    test.metrics = NULL,
    test.pred = NULL,
    cell.types = "-",
    features = "-",
    test.deconv.metrics = NULL
  ) {
    .Object@model <- model
    .Object@training.history <- training.history
    .Object@test.metrics <- test.metrics
    .Object@test.pred <- test.pred
    .Object@cell.types <- cell.types
    .Object@features <- features
    .Object@test.deconv.metrics <- test.deconv.metrics
    return(.Object)
  }
)

setMethod(
  f = "show",
  signature = "DeconvDLModel",
  definition = function(object) {
    # cat("An object of class", class(object), "\n")
    if (is.null(object@model)) {
      cat("DeconvDLModel object empty")
    } else {
      cat(paste("Trained model:", object@training.history$params$epochs,
                "epochs\n"))
      train.metrics <- lapply(object@training.history$metrics,
                              function(x) x[length(x)])
      cat("  Training metrics (last epoch):\n")
      cat(paste0("    ", names(train.metrics), ": ",
                 lapply(train.metrics, round, 4),
                 collapse = "\n"))
      cat("\n  Evaluation metrics on test data:\n")
      cat(paste0("    ", names(object@test.metrics), ": ",
                 lapply(object@test.metrics, round, 4),
                 collapse = "\n"))
      if (!is.null(test.deconv.metrics(object)) ||
          length(test.deconv.metrics(object)) > 0) {
        cat("\n  Performance evaluation over each sample: ")
        cat(paste(names(test.deconv.metrics(object)[[2]]), collapse = " "))
      }
    }
  }
)

setClassUnion("DeconvDLModelOrNULL", c("DeconvDLModel", "NULL"))


################################################################################
########################### SpatialDDLS class ##############################
################################################################################

## TODO: create a function to check if spatial.experiments contains SpatialExperiment objects

#' The SpatialDDLS Class
#'
#' The \code{\linkS4class{SpatialDDLS}}object is the core of the
#' \pkg{SpatialDDLS} package. This object stores different intermediate data
#' needed for the construction of new deconvolution models, the spatial
#' transcriptomics profiles to be deconvoluted, and the predicted cell type
#' proportions.
#'
#' This object uses other classes to store different types of data generated
#' during the workflow: \itemize{ \item
#' \code{\linkS4class{SingleCellExperiment}} class for single-cell RNA-Seq data
#' storage, using sparse matrix from the \pkg{Matrix} package
#' (\code{\linkS4class{dgCMatrix}} class) or \code{HDF5Array} class in case of
#' using HDF5 files as back-end (see below for more information). \item
#' \code{\linkS4class{SpatialExperiment}} class for spatial transcriptomics data
#' storage. \item \code{\linkS4class{ZinbModel}} class with estimated parameters
#' for the simulation of new single-cell profiles. \item
#' \code{\linkS4class{SummarizedExperiment}} class for simulated mixed
#' transcriptional profiles storage. \item \code{\linkS4class{PropCellTypes}}
#' class for compositional cell type matrices. See
#' \code{?\linkS4class{PropCellTypes}} for details. \item
#' \code{\linkS4class{DeconvDLModel}} class to store the information related to
#' deep neural network models. See \code{?\linkS4class{DeconvDLModel}} for
#' details. }
#'
#' In order to provide a way to work with large amounts of data on
#' RAM-constrained machines, we provide the possibility to use HDF5 files as
#' back-end to store count matrices of both real and simulated single-cell
#' profiles by using the \pkg{HDF5Array} and \pkg{DelayedArray} classes from the
#' homonymous packages.
#'
#' @slot single.cell.real Real single-cell data stored in a
#'   \code{SingleCellExperiment} object. The count matrix is stored either as
#'   \code{\linkS4class{dgCMatrix}} or \code{HDF5Array} objects.
#' @slot spatial.experiments List of \code{\linkS4class{SpatialExperiment}}
#'   objects to be deconvoluted.
#' @slot zinb.params \code{\linkS4class{ZinbModel}} object with estimated
#'   parameters for the simulation of new single-cell expression profiles.
#' @slot single.cell.simul Simulated single-cell expression profiles using the
#'   ZINB-WaVE model.
#' @slot prob.cell.types \code{\linkS4class{PropCellTypes}} class with cell
#'   composition matrices built for the simulation of mixed transcriptional
#'   profiles with known cell composition.
#' @slot mixed.profiles List of simulated train and test mixed transcriptional
#'   profiles. Each entry is a \code{\linkS4class{SummarizedExperiment}} object.
#'   Count matrices can be stored as \code{HDF5Array} objects using HDF5 files
#'   as back-end in case of RAM limitations.
#' @slot trained.model \code{\linkS4class{DeconvDLModel}} object with
#'   information related to the trained model. See
#'   \code{?\linkS4class{DeconvDLModel}} for more details.
#' @slot deconv.spots Deconvolution results. It consists of a list where each
#'   element corresponds to the results for each
#'   \code{\linkS4class{SpatialExperiment}} object contained in the
#'   \code{spatial.experiments} slot.
#' @slot project Name of the project.
#' @slot version Version of \pkg{SpatialDDLS} this object was built under.
#'
#' @exportClass SpatialDDLS
#' @export SpatialDDLS
#'   
SpatialDDLS <- setClass(
  Class = "SpatialDDLS",
  slots = c(
    single.cell.real = "SingleCellExperimentOrNULL",
    spatial.experiments = "ListOrNULL",
    zinb.params = "ZinbParametersModelOrNULL",
    single.cell.simul = "SingleCellExperimentOrNULL",
    prob.cell.types = "ListOrNULL",
    mixed.profiles = "ListOrNULL",
    trained.model = "DeconvDLModelOrNULL",
    deconv.spots = "ListOrNULL",
    project = "character",
    version = "package_version"
  )
)

setMethod(
  f = "initialize", 
  signature = "SpatialDDLS",
  definition = function(
    .Object,
    single.cell.real = NULL,
    spatial.experiments = NULL, 
    zinb.params = NULL,
    single.cell.simul = NULL,
    prob.cell.types = NULL,
    mixed.profiles = NULL, 
    trained.model = NULL,
    deconv.spots = NULL,
    project = "SpatialDDLSProject",
    version = packageVersion(pkg = "digitalDLSorteR")
  ) {
    .Object@single.cell.real <- single.cell.real
    .Object@spatial.experiments <- spatial.experiments
    .Object@zinb.params <- zinb.params
    .Object@single.cell.simul <- single.cell.simul
    .Object@prob.cell.types <- prob.cell.types
    .Object@mixed.profiles <- mixed.profiles
    .Object@trained.model <- trained.model
    .Object@deconv.spots <- deconv.spots
    .Object@project <- project
    .Object@version <- version
    return(.Object)
  }
)

.selectSome <- function(vec, num) {
  if (length(vec) < 6)
    namesSel <- sample(length(vec), size = 6, replace = TRUE)    
  else 
    namesSel <- sample(length(vec), size = 6)    
  return(
    paste(
      paste(vec[namesSel[1:3]], collapse = " "), "...", 
      paste(vec[namesSel[3:6]], collapse = " ")
    )
  )
}

.sceShow <- function(sce) {
  cat(" ", dim(sce)[1], "features and", dim(sce)[2], "cells\n")
  if (is.null(rownames(sce))) rownames.sce <- "---"
  else rownames.sce <- .selectSome(vec = rownames(sce), num = 6)
  if (identical(colnames(sce), character(0))) colnames.sce <- "---"
  else colnames.sce <- .selectSome(vec = colnames(sce), num = 6)
  cat("  rownames:", rownames.sce, "\n")
  cat("  colnames:", colnames.sce, "\n")
}

.bulkShow <- function(se) {
  cat("   ", dim(se)[1], "features and", dim(se)[2], "samples\n")
  if (dim(rowData(se))[2] == 0) rownames.se <- "---" 
  else rownames.se <- .selectSome(vec = rownames(rowData(se)), num = 6)
  if (identical(colnames(se), character(0))) colnames.se <- "---"
  else colnames.se <- .selectSome(vec = rownames(colData(se)), num = 6)
  cat("    rownames:", rownames.se, "\n")
  cat("    colnames:", colnames.se, "\n")
}

.finalShow <- function(se) {
  cat("   ", dim(se)[2], "features and", dim(se)[1], "samples: ")
  n.spot <- sum(grepl("Spot\\.*", rowData(se)[[1]]))
  n.sc <- abs(n.spot - dim(se)[1])
  cat(n.spot, "mixed spots and", n.sc, "single-cell profiles\n")
}

.zinbModelShow <- function(zinb.model) {
  cat(
    paste0(
      "ZinbModel object:\n",
      "  ", zinbwave::nSamples(zinb.model), " samples; ",
      "  ", zinbwave::nFeatures(zinb.model), " genes.\n",
      "  ", NCOL(zinbwave::getX_mu(zinb.model)),
      " sample-level covariate(s) (mu); ",
      "  ", NCOL(zinbwave::getX_pi(zinb.model)),
      " sample-level covariate(s) (pi);\n",
      "  ", NCOL(zinbwave::getV_mu(zinb.model)),
      " gene-level covariate(s) (mu); ",
      "  ", NCOL(zinbwave::getV_pi(zinb.model)),
      " gene-level covariate(s) (pi);\n",
      "  ", zinbwave::nFactors(zinb.model), " latent factor(s).\n"
    )
  )
}

.allSlotsNull <- function(object) {
  list.slots <- list(
    "single.cell.real", "spatial.experiments", "zinb.params", 
    "single.cell.simul", "prob.cell.types", "mixed.profiles", 
    "trained.model", "deconv.spots"
  )
  res <- all(
    unlist(
      lapply(list.slots, function(x) is.null(do.call("@", list(object, x))))
    )
  )
  if (res) return(TRUE)
  else return(FALSE)
}

setMethod(
  f = "show",
  signature = "SpatialDDLS",
  definition = function(object) {
    if (.allSlotsNull(object)) {
      cat("An empty object of class", class(object), "\n")
      opt <- options(show.error.messages = FALSE)
      on.exit(options(opt))
      stop()
    } else {
      cat("An object of class", class(object), "\n")
    }
    if (!is.null(object@single.cell.real)) {
      cat("Real single-cell profiles:\n")
      .sceShow(object@single.cell.real)
    } else {
      cat("Real single-cell profiles:\n")
      .sceShow(S4Vectors::DataFrame())
    }
    if (!is.null(object@spatial.experiments)) {
      cat("Spatial experiments:\n")
      cat(" ", length(object@spatial.experiments), "experiments\n")
      if (length(object@spatial.experiments) < 3)
        sapply(object@spatial.experiments, .sceShow)
    } else {
      cat("Spatial experiments:\n")
      .sceShow(S4Vectors::DataFrame())
    }
    if (!is.null(object@zinb.params)) {
      .zinbModelShow(object@zinb.params@zinbwave.model)
    }
    if (!is.null(object@single.cell.simul)) {
      cat("Simulated single-cell profiles:\n")
      .sceShow(object@single.cell.simul)
    }
    if (!is.null(object@prob.cell.types)) {
      cat("Cell type composition matrices:\n")
      lapply(X = c("train", "test"), FUN = function(x) {
        if (x %in% names(object@prob.cell.types)) {
          cat(show(object@prob.cell.types[[x]]), "\n")
        }
      })
    }
    if (!is.null(object@mixed.profiles)) {
      cat("Simulated mixed spots:\n")
      lapply(
        X = c("train", "test"), FUN = function(x) {
          if (x %in% names(object@mixed.profiles)) {
            cat(paste(" ", x, "spots:\n"))
            .bulkShow(object@mixed.profiles[[x]])
          }
        }
      )
    }
    if (!is.null(object@trained.model)) {
      cat(show(object@trained.model), "\n")
    }
    if (!is.null(object@deconv.spots)) {
      cat("Results (estimated cell proportions):\n")
      lapply(
        X = names(object@deconv.spots), FUN = function(x) {
          if (x %in% names(object@deconv.spots)) {
            cat(paste("  Results of", x, "mixed spots\n"))
          }
        }
      )
    }
    cat("Project:", object@project, "\n")
  }
)

.onLoad <- function(libname, pkgname) {
  if (.isConda()) {
    tryCatch(
      expr = reticulate::use_condaenv("spatialddls-env", required = TRUE),
      error = function(e) NULL
    )
  }
  Sys.setenv(TF_CPP_MIN_LOG_LEVEL = 2)
}
