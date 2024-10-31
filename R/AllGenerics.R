#' @include AllClasses.R
NULL

################################################################################
################# getters and setters for PropCellTypes class ##################
################################################################################

# prob.matrix

#' @title Get and set \code{prob.matrix} slot in a
#'   \code{\linkS4class{PropCellTypes}} object
#'
#' @docType methods
#' @name prob.matrix
#' @rdname prob.matrix
#' @aliases prob.matrix,PropCellTypes-method
#'
#' @param object \code{\linkS4class{PropCellTypes}} object.
#'
#' @export prob.matrix
#'   
setGeneric(
  name = "prob.matrix", def = function(object) standardGeneric("prob.matrix")
)
setMethod(
  f = "prob.matrix",
  signature = "PropCellTypes",
  definition = function(object) object@prob.matrix
)


#' @docType methods
#' @rdname prob.matrix
#' @aliases prob.matrix<-,PropCellTypes-method
#'
#' @param value Matrix with cell types as columns and samples as rows.
#'
#' @export prob.matrix<-
#'   
setGeneric(
  name = "prob.matrix<-", 
  def = function(object, value) standardGeneric("prob.matrix<-")
)
setMethod(
  f = "prob.matrix<-",
  signature = "PropCellTypes",
  definition = function(object, value) {
    object@prob.matrix <- value
    return(object)
  }
)

# cell.names

#' @title Get and set \code{cell.names} slot in a
#'   \code{\linkS4class{PropCellTypes}} object
#'
#' @docType methods
#' @name cell.names
#' @rdname cell.names
#' @aliases cell.names,PropCellTypes-method
#'
#' @param object \code{\linkS4class{PropCellTypes}} object.
#'
#' @export cell.names
#'   
setGeneric(
  name = "cell.names", def = function(object) standardGeneric("cell.names")
)
setMethod(
  f = "cell.names",
  signature = "PropCellTypes",
  definition = function(object) object@cell.names
)

#' @docType methods
#' @rdname cell.names
#' @aliases cell.names<-,PropCellTypes-method
#'
#' @param value Matrix containing names of the mixed transcriptional profiles to
#'   be simulated as rows and cells to be used to simulate them as columns.
#'
#' @export cell.names<-
#'   
setGeneric(
  name = "cell.names<-", 
  def = function(object, value) standardGeneric("cell.names<-")
)
setMethod(
  f = "cell.names<-",
  signature = "PropCellTypes",
  definition = function(object, value) {
    object@cell.names <- value
    return(object)
  }
)

# set.list

#' @title Get and set \code{set.list} slot in a
#'   \code{\linkS4class{PropCellTypes}} object
#'
#' @docType methods
#' @name set.list
#' @rdname set.list
#' @aliases set.list,PropCellTypes-method
#'
#' @param object \code{\linkS4class{PropCellTypes}} object.
#'
#' @export set.list
#'   
setGeneric(
  name = "set.list", def = function(object) standardGeneric("set.list")
)
setMethod(
  f = "set.list",
  signature = "PropCellTypes",
  definition = function(object) object@set.list
)

#' @docType methods
#' @rdname set.list
#' @aliases set.list<-,PropCellTypes-method
#'
#' @param value List of cells sorted by their corresponding cell type.
#'
#' @export set.list<-
#'   
setGeneric(
  name = "set.list<-", 
  def = function(object, value) standardGeneric("set.list<-")
)
setMethod(
  f = "set.list<-",
  signature = "PropCellTypes",
  definition = function(object, value) {
    object@set.list <- value
    return(object)
  }
)

# set

#' @title Get and set \code{set} slot in a \code{\linkS4class{PropCellTypes}}
#'   object
#'
#' @docType methods
#' @name set
#' @rdname set
#' @aliases set,PropCellTypes-method
#'
#' @param object \code{\linkS4class{PropCellTypes}} object.
#'
#' @export set
#'   
setGeneric(name = "set", def = function(object) standardGeneric("set"))
setMethod(
  f = "set",
  signature = "PropCellTypes",
  definition = function(object) object@set
)

#' @docType methods
#' @rdname set
#' @aliases set<-,PropCellTypes-method
#'
#' @param value A vector containing the names of cells that are present in the
#'   object.
#'
#' @export set<-
#'   
setGeneric(
  name = "set<-", def = function(object, value) standardGeneric("set<-")
)
setMethod(
  f = "set<-",
  signature = "PropCellTypes",
  definition = function(object, value) {
    object@set <- value
    return(object)
  }
)

# method

#' @title Get and set \code{method} slot in a \code{\linkS4class{PropCellTypes}}
#'   object
#'
#' @docType methods
#' @name method
#' @rdname method
#' @aliases method,PropCellTypes-method
#'
#' @param object \code{\linkS4class{PropCellTypes}} object.
#'
#' @export method
#'   
setGeneric(name = "method", def = function(object) standardGeneric("method"))
setMethod(
  f = "method",
  signature = "PropCellTypes",
  definition = function(object) object@method
)

#' @docType methods
#' @rdname method
#' @aliases method<-,PropCellTypes-method
#'
#' @param value Vector containing the method by which cell type proportions were
#'   generated.
#'
#' @export method<-
#'   
setGeneric(
  name = "method<-", def = function(object, value) standardGeneric("method<-")
)
setMethod(
  f = "method<-",
  signature = "PropCellTypes",
  definition = function(object, value) {
    object@method <- value
    return(object)
  }
)


# plots

#' @title Get and set \code{plots} slot in a \code{\linkS4class{PropCellTypes}}
#'   object
#'
#' @docType methods
#' @name plots
#' @rdname plots
#' @aliases plots,PropCellTypes-method
#'
#' @param object \code{\linkS4class{PropCellTypes}} object.
#'
#' @export plots
#'   
setGeneric(name = "plots", def = function(object) standardGeneric("plots"))
setMethod(
  f = "plots",
  signature = "PropCellTypes",
  definition = function(object) object@plots
)

#' @docType methods
#' @rdname plots
#' @aliases plots<-,PropCellTypes-method
#'
#' @param value List of lists with plots showing the distribution of cell
#'   proportions generated by each method.
#'
#' @export plots<-
#'   
setGeneric(
  name = "plots<-", def = function(object, value) standardGeneric("plots<-")
)
setMethod(
  f = "plots<-",
  signature = "PropCellTypes",
  definition = function(object, value) {
    object@plots <- value
    return(object)
  }
)

################################################################################
################# getters and setters for DeconvDLModel class ##################
################################################################################

# model

#' @title Get and set \code{model} slot in a \code{\linkS4class{DeconvDLModel}}
#'   object
#'
#' @docType methods
#' @name model
#' @rdname model
#' @aliases model,DeconvDLModel-method
#'
#' @param object \code{\linkS4class{DeconvDLModel}} object.
#'
#' @export model
#'   
setGeneric(name = "model", def = function(object) standardGeneric("model"))
setMethod(
  f = "model",
  signature = "DeconvDLModel",
  definition = function(object) object@model
)

#' @docType methods
#' @rdname model
#' @aliases model<-,DeconvDLModel-method
#'
#' @param value \code{keras.engine.sequential.Sequential} object with a trained
#'   deep neural network model.
#'
#' @export model<-
#'   
setGeneric(
  name = "model<-", def = function(object, value) standardGeneric("model<-")
)
setMethod(
  f = "model<-",
  signature = "DeconvDLModel",
  definition = function(object, value) {
    object@model <- value
    return(object)
  }
)

# training.history

#' @title Get and set \code{training.history} slot in a
#'   \code{\linkS4class{DeconvDLModel}} object
#'
#' @docType methods
#' @name training.history
#' @rdname training.history
#' @aliases training.history,DeconvDLModel-method
#'
#' @param object \code{\linkS4class{DeconvDLModel}} object.
#'
#' @export training.history
#'   
setGeneric(
  name = "training.history", 
  def = function(object) standardGeneric("training.history")
)
setMethod(
  f = "training.history",
  signature = "DeconvDLModel",
  definition = function(object) object@training.history
)

#' @docType methods
#' @rdname training.history
#' @aliases training.history<-,DeconvDLModel-method
#'
#' @param value \code{keras_training_history} object with the training history
#'   of the deep neural network model.
#'
#' @export training.history<-
#'   
setGeneric(
  name = "training.history<-", 
  def = function(object, value) standardGeneric("training.history<-")
)
setMethod(
  f = "training.history<-",
  signature = "DeconvDLModel",
  definition = function(object, value) {
    object@training.history <- value
    return(object)
  }
)

# test.metrics

#' @title Get and set \code{test.metrics} slot in a
#'   \code{\linkS4class{DeconvDLModel}} object
#'
#' @docType methods
#' @name test.metrics
#' @rdname test.metrics
#' @aliases test.metrics,DeconvDLModel-method
#'
#' @param object \code{\linkS4class{DeconvDLModel}} object.
#'
#' @export test.metrics
#'   
setGeneric(
  name = "test.metrics", def = function(object) standardGeneric("test.metrics")
)
setMethod(
  f = "test.metrics",
  signature = "DeconvDLModel",
  definition = function(object) object@test.metrics
)

#' @docType methods
#' @rdname test.metrics
#' @aliases test.metrics<-,DeconvDLModel-method
#'
#' @param value List with evaluation metrics after prediction on test data.
#'
#' @export test.metrics<-
#'   
setGeneric(
  name = "test.metrics<-", 
  def = function(object, value) standardGeneric("test.metrics<-")
)
setMethod(
  f = "test.metrics<-",
  signature = "DeconvDLModel",
  definition = function(object, value) {
    object@test.metrics <- value
    return(object)
  }
)

# test.pred

#' @title Get and set \code{test.pred} slot in a
#'   \code{\linkS4class{DeconvDLModel}} object
#'
#' @docType methods
#' @name test.pred
#' @rdname test.pred
#' @aliases test.pred,DeconvDLModel-method
#'
#' @param object \code{\linkS4class{DeconvDLModel}} object.
#'
#' @export test.pred
#'   
setGeneric(
  name = "test.pred", def = function(object) standardGeneric("test.pred")
)
setMethod(
  f = "test.pred",
  signature = "DeconvDLModel",
  definition = function(object) object@test.pred
)

#' @docType methods
#' @rdname test.pred
#' @aliases test.pred<-,DeconvDLModel-method
#'
#' @param value Matrix object with prediction results on test data.
#'
#' @export test.pred<-
#'   
setGeneric(
  name = "test.pred<-", 
  def = function(object, value) standardGeneric("test.pred<-")
)
setMethod(
  f = "test.pred<-",
  signature = "DeconvDLModel",
  definition = function(object, value) {
    object@test.pred <- value
    return(object)
  }
)

# cell.types

#' @title Get and set \code{cell.types} slot in a
#'   \code{\linkS4class{DeconvDLModel}} object
#'
#' @docType methods
#' @name cell.types
#' @rdname cell.types
#' @aliases cell.types,DeconvDLModel-method
#'
#' @param object \code{\linkS4class{DeconvDLModel}} object.
#'
#' @export cell.types
#'   
setGeneric(
  name = "cell.types", def = function(object) standardGeneric("cell.types")
)
setMethod(
  f = "cell.types",
  signature = "DeconvDLModel",
  definition = function(object) object@cell.types
)

#' @docType methods
#' @rdname cell.types
#' @aliases cell.types<-,DeconvDLModel-method
#'
#' @param value Vector with cell types considered by the deep neural network
#'   model.
#'
#' @export cell.types<-
#'   
setGeneric(
  name = "cell.types<-", 
  def = function(object, value) standardGeneric("cell.types<-")
)
setMethod(
  f = "cell.types<-",
  signature = "DeconvDLModel",
  definition = function(object, value) {
    object@cell.types <- value
    return(object)
  }
)

# features

#' @title Get and set \code{features} slot in a
#'   \code{\linkS4class{DeconvDLModel}} object
#'
#' @docType methods
#' @name features
#' @rdname features
#' @aliases features,DeconvDLModel-method
#'
#' @param object \code{\linkS4class{DeconvDLModel}} object.
#'
#' @export features
#'   
setGeneric(
  name = "features", def = function(object) standardGeneric("features")
)
setMethod(
  f = "features",
  signature = "DeconvDLModel",
  definition = function(object) object@features
)

#' @docType methods
#' @rdname features
#' @aliases features<-,DeconvDLModel-method
#'
#' @param value Vector with features (genes) considered by the deep neural
#'   network model.
#'
#' @export features<-
#'   
setGeneric(
  name = "features<-", 
  def = function(object, value) standardGeneric("features<-")
)
setMethod(
  f = "features<-",
  signature = "DeconvDLModel",
  definition = function(object, value) {
    object@features <- value
    return(object)
  }
)

# test.deconv.metrics

#' @title Get and set \code{test.deconv.metrics} slot in a
#'   \code{\linkS4class{DeconvDLModel}} object
#'
#' @docType methods
#' @name test.deconv.metrics
#' @rdname test.deconv.metrics
#' @aliases test.deconv.metrics,DeconvDLModel-method
#'
#' @param object \code{\linkS4class{DeconvDLModel}} object.
#' @param metrics Metrics to show (\code{'All'} by default)
#'
#' @export test.deconv.metrics
#'   
setGeneric(
  name = "test.deconv.metrics",
  def = function(object, metrics = "All") standardGeneric("test.deconv.metrics")
)
setMethod(
  f = "test.deconv.metrics",
  signature = "DeconvDLModel",
  definition = function(object, metrics) {
    if (metrics == "All") object@test.deconv.metrics
    else {
      if (!all(metrics %in% names(object@test.deconv.metrics)))
        stop("Metric provided is not present in DeconvDLModel object")
      return(object@test.deconv.metrics[[metrics]])
    }
  }
)

#' @docType methods
#' @rdname test.deconv.metrics
#' @aliases test.deconv.metrics<-,DeconvDLModel-method
#'
#' @param value List with evaluation metrics to assess the performance of the
#'   model on each sample of test data.
#' @export test.deconv.metrics<-
#'   
setGeneric(
  name = "test.deconv.metrics<-",
  def = function(object, metrics = "All", value) {
    standardGeneric("test.deconv.metrics<-")
  }
)
setMethod(
  f = "test.deconv.metrics<-",
  signature = "DeconvDLModel",
  definition = function(object, metrics, value) {
    if (metrics == "All") object@test.deconv.metrics <- value
    else {
      if (!all(metrics %in% names(object@test.deconv.metrics)))
        stop("Metric provided is not present in DeconvDLModel object")
      object@test.deconv.metrics[[metrics]] <- value
    }
    return(object)
  }
)

################################################################################
################## getters and setters for SpatialDDLS class ###################
################################################################################

# single.cell.real

#' @title Get and set \code{single.cell.real} slot in a
#'   \code{\link[SpatialExperiment]{SpatialExperiment}} object
#'
#' @docType methods
#' @name single.cell.real
#' @rdname single.cell.real
#' @aliases single.cell.real,SpatialDDLS-method
#'
#' @param object \code{\link[SpatialExperiment]{SpatialExperiment}} object.
#'
#' @export single.cell.real
#'   
setGeneric(
  name = "single.cell.real", 
  def = function(object) standardGeneric("single.cell.real")
)
setMethod(
  f = "single.cell.real",
  signature = "SpatialDDLS",
  definition = function(object) object@single.cell.real
)

#' @docType methods
#' @rdname single.cell.real
#' @aliases single.cell.real<-,SpatialDDLS-method
#'
#' @param value \code{\link[SingleCellExperiment]{SingleCellExperiment}} object 
#'   with real single-cell profiles.
#'
#' @export single.cell.real<-
#'   
setGeneric(
  name = "single.cell.real<-", 
  def = function(object, value) standardGeneric("single.cell.real<-")
)
setMethod(
  f = "single.cell.real<-",
  signature = "SpatialDDLS",
  definition = function(object, value) {
    object@single.cell.real <- value
    return(object)
  }
)

# spatial.experiments

#' @title Get and set \code{spatial.experiments} slot in a
#'   \code{\link[SpatialExperiment]{SpatialExperiment}} object
#'
#' @docType methods
#' @name spatial.experiments
#' @rdname spatial.experiments
#' @aliases spatial.experiments,SpatialDDLS-method
#'
#' @param object \code{\link[SpatialExperiment]{SpatialExperiment}} object.
#' @param index.st Index of the spatial transcriptomics data within the list. It
#'   can be either an position or a name if a named list was provided. If
#'   \code{NULL} (by default), all data contained in the
#'   \code{spatial.experiments} slot are returned.
#'
#' @export spatial.experiments
#'   
setGeneric(
  name = "spatial.experiments", 
  def = function(object, index.st = NULL) standardGeneric("spatial.experiments")
)
setMethod(
  f = "spatial.experiments",
  signature = "SpatialDDLS",
  definition = function(object, index.st) {
    if (is.null(index.st)) {
      return(object@spatial.experiments)
    } else {
      if (is.character(index.st) & !is.null(names(object@spatial.experiments))) {
        if (!index.st %in% names(object@spatial.experiments)) 
          stop("'index.st' provided does not exists in the spatial.experiments slot")
      } else if (is.character(index.st) & is.null(names(object@spatial.experiments))) {
        stop("spatial.experiment slot does not contain names, so `index.st` must be an integer")
      } 
    }
    if (length(index.st) > 1) {
      object@spatial.experiments[index.st]
    } else {
      object@spatial.experiments[[index.st]]
    }
  }
)

#' @docType methods
#' @rdname spatial.experiments
#' @aliases spatial.experiments<-,SpatialDDLS-method
#'
#' @param value List in which each element is a
#'   \code{\link[SpatialExperiment]{SpatialExperiment}} 
#'   object. It can be a named list.
#'
#' @export spatial.experiments<-
#'   
setGeneric(
  name = "spatial.experiments<-", 
  def = function(object, index.st = NULL, value) {
    standardGeneric("spatial.experiments<-")
  }
)
setMethod(
  f = "spatial.experiments<-",
  signature = "SpatialDDLS",
  definition = function(object, index.st, value) {
    if (is.null(index.st)) object@spatial.experiments <- value
    else {
      if (!index.st %in% names(object@spatial.experiments)) {
        warning(
          "'index.st' provided already exists in spatial.experiments slot. ", 
          "It will be overwritten"
        )
      }
      object@spatial.experiments[[index.st]] <- value
    }
    return(object)
  }
)

# zinb.params

#' @title Get and set \code{zinb.params} slot in a
#'   \code{\link[SpatialExperiment]{SpatialExperiment}} object
#'
#' @docType methods
#' @name zinb.params
#' @rdname zinb.params
#' @aliases zinb.params,SpatialDDLS-method
#'
#' @param object \code{\link[SpatialExperiment]{SpatialExperiment}} object.
#'
#' @export zinb.params
#'   
setGeneric(
  name = "zinb.params", def = function(object) standardGeneric("zinb.params")
)
setMethod(
  f = "zinb.params",
  signature = "SpatialDDLS",
  definition = function(object) object@zinb.params
)

#' @docType methods
#' @rdname zinb.params
#' @aliases zinb.params<-,SpatialDDLS-method
#'
#' @param value \code{\linkS4class{ZinbParametersModel}} object with a valid
#'   \code{\link[zinbwave]{zinbModel}} object.
#'
#' @export zinb.params<-
#'   
setGeneric(
  name = "zinb.params<-", 
  def = function(object, value) standardGeneric("zinb.params<-")
)
setMethod(
  f = "zinb.params<-",
  signature = "SpatialDDLS",
  definition = function(object, value) {
    object@zinb.params <- value
    return(object)
  }
)

# single.cell.simul

#' @title Get and set \code{single.cell.simul} slot in a
#'   \code{\link[SpatialExperiment]{SpatialExperiment}} object
#'
#' @docType methods
#' @name single.cell.simul
#' @rdname single.cell.simul
#' @aliases single.cell.simul,SpatialDDLS-method
#'
#' @param object \code{\link[SpatialExperiment]{SpatialExperiment}} object.
#'
#' @export single.cell.simul
#'   
setGeneric(
  name = "single.cell.simul", 
  def = function(object) standardGeneric("single.cell.simul")
)
setMethod(
  f = "single.cell.simul",
  signature = "SpatialDDLS",
  definition = function(object) object@single.cell.simul
)

#' @docType methods
#' @rdname single.cell.simul
#' @aliases single.cell.simul<-,SpatialDDLS-method
#'
#' @param value \code{\link[SingleCellExperiment]{SingleCellExperiment}} object 
#'   with simulated single-cell profiles.
#'
#' @export single.cell.simul<-
#'   
setGeneric(
  name = "single.cell.simul<-", 
  def = function(object, value) standardGeneric("single.cell.simul<-")
)
setMethod(
  f = "single.cell.simul<-",
  signature = "SpatialDDLS",
  definition = function(object, value) {
    object@single.cell.simul <- value
    return(object)
  }
)

# prob.cell.types

#' @title Get and set \code{prob.cell.types} slot in a
#'   \code{\link[SpatialExperiment]{SpatialExperiment}} object
#'
#' @docType methods
#' @name prob.cell.types
#' @rdname prob.cell.types
#' @aliases prob.cell.types,SpatialDDLS-method
#'
#' @param object \code{\link[SpatialExperiment]{SpatialExperiment}} object.
#' @param type.data Type of data to return. It can be \code{'both'} (default),
#'   \code{'train'}, or \code{'test'}.
#'
#' @export prob.cell.types
#'   
setGeneric(
  name = "prob.cell.types", 
  def = function(object, type.data = "both") standardGeneric("prob.cell.types")
)
setMethod(
  f = "prob.cell.types",
  signature = "SpatialDDLS",
  definition = function(object, type.data) {
    if (type.data == "train") object@prob.cell.types[["train"]]
    else if (type.data == "test") object@prob.cell.types[["test"]]
    else if (type.data == "both") object@prob.cell.types
    else stop(paste("No", type.data, "in prob.cell.types"))
  }
)

#' @docType methods
#' @rdname prob.cell.types
#' @aliases prob.cell.types<-,SpatialDDLS-method
#'
#' @param value List with two \code{\linkS4class{PropCellTypes}} objects
#'   corresponding to train and test data.
#'
#' @export prob.cell.types<-
#'   
setGeneric(
  name = "prob.cell.types<-", 
  def = function(object, type.data = "both", value) {
    standardGeneric("prob.cell.types<-")
  }
)
setMethod(
  f = "prob.cell.types<-",
  signature = "SpatialDDLS",
  definition = function(object, type.data, value) {
    if (type.data == "train") object@prob.cell.types[["train"]] <- value
    else if (type.data == "test") object@prob.cell.types[["test"]] <- value
    else if (type.data == "both") object@prob.cell.types <- value
    else stop(paste("No", type.data, "in prob.cell.types slot"))
    return(object)
  }
)

# mixed.profiles

#' Get and set \code{mixed.profiles} slot in a 
#' \code{\link[SpatialExperiment]{SpatialExperiment}} object
#'
#' @docType methods
#' @name mixed.profiles
#' @rdname mixed.profiles
#' @aliases mixed.profiles,SpatialDDLS-method
#'
#' @param object \code{\link[SpatialExperiment]{SpatialExperiment}} object.
#' @param type.data Type of data to return. It can be \code{'both'} (default),
#'   \code{'train'}, or \code{'test'}.
#'
#' @export mixed.profiles
#'   
setGeneric(
  name = "mixed.profiles", 
  def = function(object, type.data = "both") standardGeneric("mixed.profiles")
)
setMethod(
  f = "mixed.profiles",
  signature = "SpatialDDLS",
  definition = function(object, type.data) {
    if (type.data == "train") object@mixed.profiles[["train"]]
    else if (type.data == "test") object@mixed.profiles[["test"]]
    else if (type.data == "both") object@mixed.profiles
    else stop(paste("No", type.data, "in mixed.profiles slot"))
  }
)

#' @docType methods
#' @rdname mixed.profiles
#' @aliases mixed.profiles<-,SpatialDDLS-method
#'
#' @param value List with two 
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} objects, train and
#'   test, each one containing simulated mixed transcriptional profiles.
#'
#' @export mixed.profiles<-
#'   
setGeneric(
  name = "mixed.profiles<-", 
  def = function(object, type.data = "both", value) {
    standardGeneric("mixed.profiles<-")
  }
)
setMethod(
  f = "mixed.profiles<-",
  signature = "SpatialDDLS",
  definition = function(object, type.data, value) {
    if (type.data == "train") object@mixed.profiles[["train"]] <- value
    else if (type.data == "test") object@mixed.profiles[["test"]] <- value
    else if (type.data == "both") object@mixed.profiles <- value
    else stop(paste("No", type.data, "in mixed.profiles slot"))
    return(object)
  }
)

# trained.model

#' @title Get and set \code{trained.model} slot in a
#'   \code{\link[SpatialExperiment]{SpatialExperiment}} object
#'
#' @docType methods
#' @name trained.model
#' @rdname trained.model
#' @aliases trained.model,SpatialDDLS-method
#'
#' @param object \code{\link[SpatialExperiment]{SpatialExperiment}} object.
#'
#' @export trained.model
#'   
setGeneric(
  name = "trained.model", 
  def = function(object) standardGeneric("trained.model")
)
setMethod(
  f = "trained.model",
  signature = "SpatialDDLS",
  definition = function(object) object@trained.model
)

#' @docType methods
#' @rdname trained.model
#' @aliases trained.model<-,SpatialDDLS-method
#'
#' @param value \code{\linkS4class{DeconvDLModel}} object.
#'
#' @export trained.model<-
#'   
setGeneric(
  name = "trained.model<-", 
  def = function(object, value) standardGeneric("trained.model<-")
)
setMethod(
  f = "trained.model<-",
  signature = "SpatialDDLS",
  definition = function(object, value) {
    object@trained.model <- value
    return(object)
  }
)

# deconv.spots

#' @title Get and set \code{deconv.spots} slot in a
#'   \code{\link[SpatialExperiment]{SpatialExperiment}} object
#'
#' @docType methods
#' @name deconv.spots
#' @rdname deconv.spots
#' @aliases deconv.spots,SpatialDDLS-method
#'
#' @param object \code{\link[SpatialExperiment]{SpatialExperiment}} object.
#' @param index.st Name or index of predicted cell proportions (same as for the
#'   \code{spatial.experiments} slot). If \code{NULL} (by default), all
#'   results are returned.
#'
#' @export deconv.spots
#'   
setGeneric(
  name = "deconv.spots", 
  def = function(object, index.st = NULL) standardGeneric("deconv.spots")
)
setMethod(
  f = "deconv.spots",
  signature = "SpatialDDLS",
  definition = function(object, index.st) {
    if (is.null(index.st)) {
      return(object@deconv.spots)
    } else {
      if (is.character(index.st) & !is.null(names(object@deconv.spots))) {
        if (!index.st %in% names(object@deconv.spots)) 
          stop("'index.st' provided does not exists in `deconv.spots` slot")
      } else if (is.character(index.st) & is.null(names(object@deconv.spots))) {
        stop("`deconv.spots` slot does not contain names, so `index.st` must be an integer")
      } 
    }
    if (length(index.st) > 1) {
      object@deconv.spots[index.st]
    } else {
      object@deconv.spots[[index.st]]
    }
  }
)

#' @docType methods
#' @rdname deconv.spots
#' @aliases deconv.spots<-,SpatialDDLS-method
#'
#' @param value List of predicted cell type proportions for the experiments
#'   stored in the \code{spatial.experiments} slot.
#'
#' @export deconv.spots<-
#'   
setGeneric(
  name = "deconv.spots<-", 
  def = function(object, index.st = NULL, value) {
    standardGeneric("deconv.spots<-")
  }
)
setMethod(
  f = "deconv.spots<-",
  signature = "SpatialDDLS",
  definition = function(object, index.st, value) {
    if (is.null(index.st)) {
      object@deconv.spots <- value  
    } else {
      object@deconv.spots[[index.st]] <- value
    }
    return(object)
  }
)

# project

#' @title Get and set \code{project} slot in a 
#'   \code{\link[SpatialExperiment]{SpatialExperiment}} object
#'
#' @docType methods
#' @name project
#' @rdname project
#' @aliases project,SpatialDDLS-method
#'
#' @param object \code{\link[SpatialExperiment]{SpatialExperiment}} object.
#'
#' @export project
#'   
setGeneric(
  name = "project", def = function(object) standardGeneric("project")
)
setMethod(
  f = "project",
  signature = "SpatialDDLS",
  definition = function(object) object@project
)

#' @docType methods
#' @rdname project
#' @aliases project<-,SpatialDDLS-method
#'
#' @param value Character indicating the name of the project.
#'
#' @export project<-
#'   
setGeneric(
  name = "project<-", def = function(object, value) standardGeneric("project<-")
)
setMethod(
  f = "project<-",
  signature = "SpatialDDLS",
  definition = function(object, value) {
    object@project <- value
    return(object)
  }
)

################################################################################
############## getters and setters for ZinbParametersModel class ###############
################################################################################

# zinbwave.model

#' @title Get and set \code{zinbwave.model} slot in a
#'   \code{\linkS4class{ZinbParametersModel}} object
#'
#' @docType methods
#' @name zinbwave.model
#' @rdname zinbwave.model
#' @aliases zinbwave.model,ZinbParametersModel-method
#'
#' @param object \code{\linkS4class{ZinbParametersModel}} object.
#'
#' @export zinbwave.model
#'   
setGeneric(
  name = "zinbwave.model", 
  def = function(object) standardGeneric("zinbwave.model")
)
setMethod(
  f = "zinbwave.model",
  signature = "ZinbParametersModel",
  definition = function(object) object@zinbwave.model
)

#' @docType methods
#' @rdname zinbwave.model
#' @aliases zinbwave.model<-,ZinbParametersModel-method
#'
#' @param value \code{\link[zinbwave]{zinbModel}} object with the 
#'   estimated parameters to simulate new single-cell profiles.
#'
#' @export zinbwave.model<-
#'   
setGeneric(
  name = "zinbwave.model<-", 
  def = function(object, value) standardGeneric("zinbwave.model<-")
)
setMethod(
  f = "zinbwave.model<-",
  signature = "ZinbParametersModel",
  definition = function(object, value) {
    object@zinbwave.model <- value
    return(object)
  }
)


#' Save \code{\link[SpatialExperiment]{SpatialExperiment}} objects as RDS files
#'
#' Save \code{\link[SpatialExperiment]{SpatialExperiment}} and 
#' \code{\linkS4class{DeconvDLModel}} objects as RDS files. \pkg{keras} models 
#' cannot be stored natively as R objects (e.g. RData or RDS files). By saving 
#' the architecture as a JSON-like character object and the weights as a list, 
#' it is possible to retrieve a functional model and make new predictions. If 
#' the \code{trained.model} slot is empty, the function will behave as usual. 
#' \strong{Note:} with this option, the state of optimizer is not saved, only 
#' model's architecture and weights. It is possible to save the entire model as 
#' an HDF5 file with the \code{\link{saveTrainedModelAsH5}} function and load it
#' into a \code{\link[SpatialExperiment]{SpatialExperiment}} object with the
#' \code{\link{loadTrainedModelFromH5}} function. See documentation for details.
#'
#' @docType methods
#' @name saveRDS
#' @rdname saveRDS
#' @aliases saveRDS,saveRDS-method
#'
#' @param object \code{\link[SpatialExperiment]{SpatialExperiment}} or
#'   \code{\linkS4class{DeconvDLModel}} object to be saved
#' @param file File path where the object will be saved
#' @inheritParams base::saveRDS
#'
#' @return No return value, saves a 
#'   \code{\link[SpatialExperiment]{SpatialExperiment}} object as an RDS file on
#'   disk.
#'
#' @export
#'
#' @seealso \code{\link[SpatialExperiment]{SpatialExperiment}} 
#'   \code{\link{saveTrainedModelAsH5}}
#'   
setGeneric(
  name = "saveRDS", 
  def = function(
    object,
    file,
    ascii = FALSE,
    version = NULL,
    compress = TRUE,
    refhook = NULL
  ) {
    standardGeneric("saveRDS")
  }
)

#' @export
#'
#' @rdname saveRDS
setMethod(
  f = "saveRDS", 
  signature = "DeconvDLModel", 
  definition = function(
    object,
    file,
    ascii,
    version,
    compress,
    refhook
  ) {
    if ("keras.engine.sequential.Sequential" %in% class(model(object))) {
      object <- .saveModelToJSON(object)
      base::saveRDS(
        object = object,
        file = file,
        ascii = ascii,
        version = version,
        compress = compress,
        refhook = refhook
      )
    } else if (is(model(object), "list")) {
      base::saveRDS(
        object = object,
        file = file,
        ascii = ascii,
        version = version,
        compress = compress,
        refhook = refhook
      )
    } else {
      stop("No valid DeconvDLModel object")
    }
  }
)

#' @export
#'
#' @rdname saveRDS
setMethod(
  f = "saveRDS", 
  signature = "SpatialDDLS", 
  definition = function(
    object,
    file,
    ascii,
    version,
    compress,
    refhook
  ) {
    if (!is.null(trained.model(object))) {
      if ("keras.engine.sequential.Sequential" %in% 
          class(trained.model(object)@model)) {
        model.object <- .saveModelToJSON(trained.model(object))
        trained.model(object) <- model.object
      }
    }
    base::saveRDS(
      object = object,
      file = file,
      ascii = ascii,
      version = version,
      compress = compress,
      refhook = refhook
    )
  }
)
