#' @include AllClasses.R
NULL

################################################################################
############## getters and setters for PropCellTypes class ###############
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
#' @param value Matrix with cell types as columns and samples as
#'   rows.
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
#' @param value Matrix containing the name of the pseudo-bulk samples to be
#'   simulated as rows and the cells to be used to simulate them as columns
#'   (\code{n.cell} argument)
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
#' @param value List of cells sorted according to the cell type to which they
#'   belong.
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

#' @title Get and set \code{set} slot in a
#'   \code{\linkS4class{PropCellTypes}} object
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
#' @param value Vector with names of cells present in the object.
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

#' @title Get and set \code{method} slot in a
#'   \code{\linkS4class{PropCellTypes}} object
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
#' @param value Vector with names of cells present in the object.
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

#' @title Get and set \code{plots} slot in a
#'   \code{\linkS4class{PropCellTypes}} object
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
#' @param value List of lists with plots showing the distribution of the cell
#'   proportions generated by each method during the process.
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
############## getters and setters for DeconvDLModel class ################
################################################################################

# model

#' @title Get and set \code{model} slot in a
#'   \code{\linkS4class{DeconvDLModel}} object
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
#' @param value \code{keras.engine.sequential.Sequential} object with a
#' trained Deep Neural Network model.
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
#'   of the Deep Neural Network model
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
#' @param value List object with the resulting metrics after prediction
#'   on test data with the Deep Neural Network model.
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
#' @param value Matrix object with the prediction results on test data.
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
#' @param value Vector with cell types considered by the Deep Neural Network
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
#' @param value Vector with features (genes) considered by the Deep Neural
#'   Network model.
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
#' @param value List with evaluation metrics used to assess the
#'   performance of the model on each sample of test data.
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
################ getters and setters for SpatialDDLS class #################
################################################################################

# single.cell.real

#' @title Get and set \code{single.cell.real} slot in a
#'   \code{\linkS4class{SpatialDDLS}} object
#'
#' @docType methods
#' @name single.cell.real
#' @rdname single.cell.real
#' @aliases single.cell.real,SpatialDDLS-method
#' 
#' @param object \code{\linkS4class{SpatialDDLS}} object.
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
#' @param value \code{\linkS4class{SingleCellExperiment}} object with real
#'   single-cell profiles.
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
#'   \code{\linkS4class{SpatialDDLS}} object
#'
#' @docType methods
#' @name spatial.experiments
#' @rdname spatial.experiments
#' @aliases spatial.experiments,SpatialDDLS-method
#' 
#' @param object \code{\linkS4class{SpatialDDLS}} object.
#' @param name.data Name of the ST data. If \code{NULL} (by default), all
#'   data contained in the \code{spatial.experiments} slot are returned.
#'
#' @export spatial.experiments
#'   
setGeneric(
  name = "spatial.experiments", 
  def = function(object, name.data = NULL) standardGeneric("spatial.experiments")
)
setMethod(
  f = "spatial.experiments",
  signature = "SpatialDDLS",
  definition = function(object, name.data) {
    if (is.null(name.data)) object@spatial.experiments
    else {
      if (!name.data %in% names(object@spatial.experiments)) {
        stop("'name.data' provided does not exists in spatial.experiments slot")
      }
      return(object@spatial.experiments[[name.data]])
    }
  }
)

#' @docType methods
#' @rdname spatial.experiments
#' @aliases spatial.experiments<-,SpatialDDLS-method
#' 
#' @param value List whose names are the reference of the stored data.
#' 
#' @export spatial.experiments<-
#'
setGeneric(
  name = "spatial.experiments<-", 
  def = function(object, name.data = NULL, value) {
    standardGeneric("spatial.experiments<-")
  }
)
setMethod(
  f = "spatial.experiments<-",
  signature = "SpatialDDLS",
  definition = function(object, name.data, value) {
    if (is.null(name.data)) object@spatial.experiments <- value
    else {
      if (!name.data %in% names(object@spatial.experiments)) {
        warning(
          "'name.data' provided already exists in spatial.experiments slot. ", 
          "It will be overwritten"
        )
      }
      object@spatial.experiments[[name.data]] <- value
    }
    return(object)
  }
)

# zinb.params

#' @title Get and set \code{zinb.params} slot in a
#'   \code{\linkS4class{SpatialDDLS}} object
#'
#' @docType methods
#' @name zinb.params
#' @rdname zinb.params
#' @aliases zinb.params,SpatialDDLS-method
#' 
#' @param object \code{\linkS4class{SpatialDDLS}} object.
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
#'   \code{\linkS4class{ZinbModel}} object.
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
#'   \code{\linkS4class{SpatialDDLS}} object
#'
#' @docType methods
#' @name single.cell.simul
#' @rdname single.cell.simul
#' @aliases single.cell.simul,SpatialDDLS-method
#' 
#' @param object \code{\linkS4class{SpatialDDLS}} object.
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
#' @param value \code{\linkS4class{SingleCellExperiment}} object with simulated
#'   single-cell profiles.
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
#'   \code{\linkS4class{SpatialDDLS}} object
#'
#' @docType methods
#' @name prob.cell.types
#' @rdname prob.cell.types
#' @aliases prob.cell.types,SpatialDDLS-method
#' 
#' @param object \code{\linkS4class{SpatialDDLS}} object.
#' @param type.data Element of the list. Can be \code{'train'}, \code{'test'} or
#'   \code{'both'} (the last by default).
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
#' @param value List with two elements, train and test, each one with a
#'   \code{\linkS4class{PropCellTypes}} object.
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

# mixed.spot.profiles

#' Get and set \code{mixed.spot.profiles} slot in a
#'   \code{\linkS4class{SpatialDDLS}} object
#'
#' @docType methods
#' @name mixed.spot.profiles
#' @rdname mixed.spot.profiles
#' @aliases mixed.spot.profiles,SpatialDDLS-method
#' 
#' @param object \code{\linkS4class{SpatialDDLS}} object.
#' @param type.data Element of the list. Can be \code{'train'}, \code{'test'} or
#'   \code{'both'} (the last by default).
#'
#' @export mixed.spot.profiles
#'   
setGeneric(
  name = "mixed.spot.profiles", 
  def = function(object, type.data = "both") standardGeneric("mixed.spot.profiles")
)
setMethod(
  f = "mixed.spot.profiles",
  signature = "SpatialDDLS",
  definition = function(object, type.data) {
    if (type.data == "train") object@mixed.spot.profiles[["train"]]
    else if (type.data == "test") object@mixed.spot.profiles[["test"]]
    else if (type.data == "both") object@mixed.spot.profiles
    else stop(paste("No", type.data, "in mixed.spot.profiles slot"))
  }
)

#' @docType methods
#' @rdname mixed.spot.profiles
#' @aliases mixed.spot.profiles<-,SpatialDDLS-method
#' 
#' @param value List with two elements, train and test, each one being
#'   a \code{\linkS4class{SummarizedExperiment}} object with simulated bulk
#'   RNA-Seq samples.
#'
#' @export mixed.spot.profiles<-
#'   
setGeneric(
  name = "mixed.spot.profiles<-", 
  def = function(object, type.data = "both", value) {
    standardGeneric("mixed.spot.profiles<-")
  }
)
setMethod(
  f = "mixed.spot.profiles<-",
  signature = "SpatialDDLS",
  definition = function(object, type.data, value) {
    if (type.data == "train") object@mixed.spot.profiles[["train"]] <- value
    else if (type.data == "test") object@mixed.spot.profiles[["test"]] <- value
    else if (type.data == "both") object@mixed.spot.profiles <- value
    else stop(paste("No", type.data, "in mixed.spot.profiles slot"))
    return(object)
  }
)

# trained.model

#' @title Get and set \code{trained.model} slot in a
#'   \code{\linkS4class{SpatialDDLS}} object
#'
#' @docType methods
#' @name trained.model
#' @rdname trained.model
#' @aliases trained.model,SpatialDDLS-method
#'
#' @param object \code{\linkS4class{SpatialDDLS}} object.
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
#'   \code{\linkS4class{SpatialDDLS}} object
#'
#' @docType methods
#' @name deconv.spots
#' @rdname deconv.spots
#' @aliases deconv.spots,SpatialDDLS-method
#' 
#' @param object \code{\linkS4class{SpatialDDLS}} object.
#' @param name.data Name of the data. If \code{NULL} (by default), all
#'   results contained in the \code{deconv.spots} slot are returned.
#'
#' @export deconv.spots
#'   
setGeneric(
  name = "deconv.spots", 
  def = function(object, name.data = NULL) standardGeneric("deconv.spots")
)
setMethod(
  f = "deconv.spots",
  signature = "SpatialDDLS",
  definition = function(object, name.data) {
    if (is.null(name.data)) object@deconv.spots
    else {
      if (!name.data %in% names(object@deconv.spots)) {
        stop("'name.data' provided does not exists in deconv.spots slot")
      }
      return(object@deconv.spots[[name.data]])
    }
  }
)

#' @docType methods
#' @rdname deconv.spots
#' @aliases deconv.spots<-,SpatialDDLS-method
#'
#' @param value List whose names are the reference of the stored results.
#'
#' @export deconv.spots<-
#'   
setGeneric(
  name = "deconv.spots<-", 
  def = function(object, name.data = NULL, value) {
    standardGeneric("deconv.spots<-")
  }
)
setMethod(
  f = "deconv.spots<-",
  signature = "SpatialDDLS",
  definition = function(object, name.data, value) {
    if (is.null(name.data)) {
      object@deconv.spots <- value  
    } else {
      object@deconv.spots[[name.data]] <- value
    }
    return(object)
  }
)

# project

#' @title Get and set \code{project} slot in a
#'   \code{\linkS4class{SpatialDDLS}} object
#'
#' @docType methods
#' @name project
#' @rdname project
#' @aliases project,SpatialDDLS-method
#' 
#' @param object \code{\linkS4class{SpatialDDLS}} object.
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
#' @param value \code{\linkS4class{ZinbModel}} object with the estimated
#'   parameters.
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


#' Save \code{\linkS4class{SpatialDDLS}} objects as RDS files
#'
#' Save \code{\linkS4class{SpatialDDLS}} and
#' \code{\linkS4class{DeconvDLModel}} objects as RDS files. \pkg{keras}
#' models cannot be stored natively as R objects (e.g. RData or RDS files). By
#' saving the structure as a JSON-like character object and the weights as a
#' list, it is possible to retrieve the model and make predictions. If the
#' \code{trained.model} slot is empty, the function will behave as usual.
#' \strong{Note:} with this option, the state of optimizer is not saved, only
#' the architecture and weights. It is possible to save the entire model as an
#' HDF5 file with the \code{\link{saveTrainedModelAsH5}} function and to load it
#' into a \code{\linkS4class{SpatialDDLS}} object with the
#' \code{\link{loadTrainedModelFromH5}} function. See documentation for details.
#'
#' @docType methods
#' @name saveRDS
#' @rdname saveRDS
#' @aliases saveRDS,saveRDS-method
#'
#' @param object \code{\linkS4class{SpatialDDLS}} or
#'   \code{\linkS4class{DeconvDLModel}} object to be saved
#' @param file File path where the object will be saved
#' @inheritParams base::saveRDS
#'
#' @return No return value, saves a \code{\linkS4class{SpatialDDLS}} object
#'   as an RDS file on disk.
#'
#' @export
#'
#' @seealso \code{\linkS4class{SpatialDDLS}}
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

#' Bar plot of deconvoluted cell type proportions in bulk RNA-Seq samples
#'
#' Bar plot of deconvoluted cell type proportions in bulk RNA-Seq samples.
#'
#' @param data \code{\linkS4class{SpatialDDLS}} object with
#'   \code{deconv.spots} slot or a data frame/matrix with cell types as
#'   columns and samples as rows.
#' @param colors Vector of colors to be used.
#' @param simplify Type of simplification performed during deconvolution. Can be
#'   \code{simpli.set} or \code{simpli.maj} (\code{NULL} by default). It is only
#'   for \code{\linkS4class{SpatialDDLS}} objects.
#' @param color.line Color of the border bars.
#' @param x.label Label of x-axis.
#' @param rm.x.text Logical value indicating whether to remove x-axis ticks
#'   (name of samples).
#' @param title Title of the plot.
#' @param legend.title Title of the legend plot.
#' @param angle Angle of text ticks.
#' @param name.data If a \code{\linkS4class{SpatialDDLS}} is given, name of
#'   the element that stores the results in the \code{deconv.spots} slot.
#' @param theme \pkg{ggplot2} theme.
#' @param ... Other arguments for specific methods.
#'
#' @return A ggplot object with the provided cell proportions represented as a
#'   bar plot.
#'
#' @export
#'
#' @examples
#' # matrix of simulated proportions (same estructure as deconvolution results)
#' deconvResults <- gtools::rdirichlet(n = 20, alpha = c(1, 1, 1, 0.5, 0.1))
#' colnames(deconvResults) <- paste("CellType", seq(ncol(deconvResults)))
#' rownames(deconvResults) <- paste("BulkSample", seq(nrow(deconvResults)))
#' barPlotCellTypes(deconvResults)
#'
#' # Using a SpatialDDLS object
#' DDLS <- SpatialDDLS(deconv.spots = list(Example = deconvResults))
#' barPlotCellTypes(DDLS)
#' 
#' @rdname barPlotCellTypes
#'
#' @seealso \code{\link{deconvSpatialDDLS}}
#'   \code{\link{deconvSpatialDDLSObj}}
#'   
setGeneric(
  name = "barPlotCellTypes", 
  def = function(
    data,
    colors = NULL,
    simplify = NULL,
    color.line = NA,
    x.label = "Bulk samples",
    rm.x.text = FALSE,
    title = "Results of deconvolution",
    legend.title = "Cell types",
    angle = 90,
    theme = NULL, 
    ...
  ) {
    standardGeneric("barPlotCellTypes")
  }
)

#' @export
#'
#' @rdname barPlotCellTypes
setMethod(
  f = "barPlotCellTypes",
  signature = signature(data = "SpatialDDLS"),
  definition = function(
    data,
    colors = NULL,
    simplify = NULL,
    color.line = NA,
    x.label = "Bulk samples",
    rm.x.text = FALSE,
    title = "Results of deconvolution",
    legend.title = "Cell types",
    angle = 90,
    theme = NULL,
    name.data = NULL
  ) {
    if (is.null(deconv.spots(data))) {
      stop("There are no results in SpatialDDLS object. Please see ?deconvDigitalDLSorterObj")
    } else if (is.null(name.data)) {
      message("'name.data' not provided. By default, first results are used")
      name.data <- 1
    } else if (!any(name.data %in% names(deconv.spots(data))) &&
               !any(name.data %in% seq_along(deconv.spots(data)))) {
      stop("Provided 'name.data' does not exist")
    }
    if (!is.null(simplify) && !is.na(simplify)) {
      if (!is(deconv.spots(data)[[name.data]], "list")) {
        stop("No simplified results available")
      } else {
        if (simplify != "simpli.set" && simplify != "simpli.majority") {
          stop("simplify argument must be one of the following options: ",
               "'simpli.set' or 'simpli.majority'")
        } else if (!any(simplify == names(deconv.spots(data)[[name.data]]))) {
          stop(paste(simplify, "data is not present in DeconvDLModel object"))
        }
        res <- deconv.spots(data)[[name.data]][[simplify]]
      }
    } else {
      if (is(deconv.spots(data)[[name.data]], "list")) {
        res <- deconv.spots(data)[[name.data]][[1]]
      } else {
        res <- deconv.spots(data)[[name.data]]
      }
    }
    if (is.null(colnames(res))) {
      stop("'data' must have colnames (corresponding cell types). Please run deconvDigitalDLSorterObj")
    }
    return(
      .barPlot(
        data = res,
        colors = colors,
        color.line = color.line,
        x.label = x.label,
        rm.x.text = rm.x.text,
        title = title,
        legend.title = legend.title,
        angle = angle,
        theme = theme
      )
    )
  }
)

#' @export
#'
#' @rdname barPlotCellTypes
setMethod(
  f = "barPlotCellTypes", 
  signature = signature(data = "ANY"),
  definition = function(
    data,
    colors,
    color.line = NA,
    x.label = "Bulk samples",
    rm.x.text = FALSE,
    title = "Results of deconvolution",
    legend.title = "Cell types",
    angle = 90,
    theme = NULL
  ) {
    if (is.null(colnames(data))) {
      stop("'data' must have colnames (corresponding cell types). Please run deconvDigitalDLSorter")
    }
    plot <- .barPlot(
      data = data,
      colors = colors,
      color.line = color.line,
      x.label = x.label,
      rm.x.text = rm.x.text,
      title = title,
      legend.title = legend.title,
      angle = angle,
      theme = theme
    )
    return(plot)
  }
)


