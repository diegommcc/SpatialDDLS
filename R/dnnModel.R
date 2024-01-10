#' @importFrom dplyr %>%
#' @importFrom stats dist
NULL

################################################################################
################### Train and evaluate deconvolution model #####################
################################################################################

#' Train deconvolution model for spatial transcriptomics data
#'
#' Train a deep neural network model using training data from the
#' \code{\linkS4class{SpatialDDLS}} object. This model will be used to 
#' deconvolute spatial transcriptomics data from the same biological context as 
#' the single-cell RNA-seq data used to train it. In addition, the trained 
#' model is evaluated using test data, and prediction results are obtained to 
#' determine its performance (see \code{?\link{calculateEvalMetrics}}).
#'
#' \strong{Simulation of mixed transcriptional profiles 'on the fly'}
#'
#' \code{trainDeconvModel} can avoid storing simulated mixed spot profiles by
#' using the \code{on.the.fly} argument. This functionality aims at reducing the
#' the \code{simMixedProfiles} function's memory usage: simulated profiles are
#' built in each batch during training/evaluation.
#'
#' \strong{Neural network architecture}
#'
#' It is possible to change the model's architecture: number of hidden layers,
#' number of neurons for each hidden layer, dropout rate, activation function,
#' and loss function. For more customized models, it is possible to provide a
#' pre-built model through the \code{custom.model} argument (a
#' \code{keras.engine.sequential.Sequential} object) where it is necessary that
#' the number of input neurons is equal to the number of considered
#' features/genes, and the number of output neurons is equal to the number of
#' considered cell types.
#'
#' @param object \code{\linkS4class{SpatialDDLS}} object with
#'   \code{single.cell.real}/\code{single.cell.simul}, \code{prob.cell.types},
#'   and \code{mixed.profiles} slots (the last only if \code{on.the.fly =
#'   FALSE}).
#' @param type.data.train Type of profiles to be used for training. It can be
#'   \code{'both'}, \code{'single-cell'} or \code{'mixed'} (\code{'mixed'} by
#'   default).
#' @param type.data.test Type of profiles to be used for evaluation. It can be
#'   \code{'both'}, \code{'single-cell'} or \code{'mixed'} (\code{'mixed'} by
#'   default).
#' @param batch.size Number of samples per gradient update (64 by default).
#' @param num.epochs Number of epochs to train the model (60 by default).
#' @param num.hidden.layers Number of hidden layers of the neural network (2 by
#'   default). This number must be equal to the length of \code{num.units}
#'   argument.
#' @param num.units Vector indicating the number of neurons per hidden layer
#'   (\code{c(200, 200)} by default). The length of this vector must be equal to
#'   the \code{num.hidden.layers} argument.
#' @param activation.fun Activation function (\code{'relu'} by default). See
#'   the
#'   \href{https://tensorflow.rstudio.com/reference/keras/activation_relu.html}{keras
#'    documentation} to know available activation functions.
#' @param dropout.rate Float between 0 and 1 indicating the fraction of
#'   input neurons to be dropped in layer dropouts (0.25 by default). By
#'   default, \pkg{SpatialDDLS} implements 1 dropout layer per hidden layer.
#' @param loss Character indicating loss function selected for model training
#'   (\code{'kullback_leibler_divergence'} by default). See the
#'   \href{https://tensorflow.rstudio.com/reference/keras/loss-functions.html}{keras
#'    documentation} to know available loss functions.
#' @param metrics Vector of metrics used to assess model performance during
#'   training and evaluation (\code{c("accuracy", "mean_absolute_error",
#'   "categorical_accuracy")} by default). See the
#'   \href{https://tensorflow.rstudio.com/reference/keras/metric_binary_accuracy.html}{keras
#'    documentation} to know available performance metrics.
#' @param normalize Whether to normalize data using logCPM (\code{TRUE} by 
#'   default). This parameter is only considered when the method used to 
#'   simulate mixed transcriptional profiles (\code{simMixedProfiles} 
#'   function) was \code{"AddRawCount"}. Otherwise, data were already 
#'   normalized.
#' @param scaling How to scale data before training. It can be:
#'   \code{"standardize"} (values are centered around the mean with a unit
#'   standard deviation), \code{"rescale"} (values are shifted and rescaled so
#'   that they end up ranging between 0 and 1) or \code{"none"} (no
#'   scaling is performed). \code{"standardize"} by default. 
#' @param norm.batch.layers Whether to include batch normalization layers
#'   between each hidden dense layer (\code{TRUE} by default).
#' @param custom.model It allows to use a custom neural network architecture. It 
#'   must be a \code{keras.engine.sequential.Sequential} object in which the 
#'   number of input neurons is equal to the number of considered 
#'   features/genes, and the number of output neurons is equal to the number of 
#'   cell types considered (\code{NULL} by default). If provided, the arguments 
#'   related to the neural network architecture will be ignored.
#' @param shuffle Boolean indicating whether data will be shuffled (\code{TRUE}
#'   by default).
#' @param sc.downsampling It is only used if \code{type.data.train} is equal to
#'   \code{'both'} or \code{'single-cell'}. It allows to set a maximum number of
#'   single-cell profiles of a specific cell type for training to avoid
#'   an unbalanced representation of classes (\code{NULL} by default).
#' @param use.generator Boolean indicating whether to use generators during
#'   training and test. Generators are automatically used when \code{on.the.fly
#'   = TRUE} or HDF5 files are used, but it can be activated by the user on
#'   demand (\code{FALSE} by default).
#' @param on.the.fly Boolean indicating whether simulated data will be generated
#'   'on the fly' during training (\code{FALSE} by default).
#' @param agg.function If \code{on.the.fly == TRUE}, function used to build 
#'   mixed transcriptional profiles. It may be: 
#'   \itemize{ \item \code{"AddRawCount"} (by default): single-cell
#'   profiles (raw counts) are added up across cells. Then, log-CPMs are
#'   calculated. \item \code{"MeanCPM"}: single-cell profiles (raw counts) are
#'   transformed into logCPM and cross-cell averages are calculated. 
#'   \item \code{"AddCPM"}: single-cell profiles (raw counts) are transformed 
#'   into CPMs and are added up across cells. Then, log-CPMs are calculated.}
#' @param threads Number of threads used during simulation of mixed
#'   transcriptional profiles if \code{on.the.fly = TRUE} (1 by default).
#' @param view.metrics.plot Boolean indicating whether to show plots of loss and
#'   evaluation metrics during training (\code{TRUE} by default). \pkg{keras}
#'   for R allows to see model progression during training if you are working in
#'   RStudio.
#' @param verbose Boolean indicating whether to display model progression during
#'   training and model architecture information (\code{TRUE} by default).
#'
#' @return A \code{\linkS4class{SpatialDDLS}} object with \code{trained.model}
#'   slot containing a \code{\linkS4class{DeconvDLModel}} object. For more
#'   information about the structure of this class, see
#'   \code{?\linkS4class{DeconvDLModel}}.
#'
#' @export
#'
#' @seealso \code{\link{plotTrainingHistory}} \code{\link{deconvSpatialDDLS}}
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(
#'     counts = matrix(
#'       rpois(30, lambda = 5), nrow = 15, ncol = 10,
#'       dimnames = list(paste0("Gene", seq(15)), paste0("RHC", seq(10)))
#'     )
#'   ),
#'   colData = data.frame(
#'     Cell_ID = paste0("RHC", seq(10)),
#'     Cell_Type = sample(x = paste0("CellType", seq(2)), size = 10,
#'                        replace = TRUE)
#'   ),
#'   rowData = data.frame(
#'     Gene_ID = paste0("Gene", seq(15))
#'   )
#' )
#' SDDLS <- createSpatialDDLSobject(
#'   sc.data = sce,
#'   sc.cell.ID.column = "Cell_ID",
#'   sc.gene.ID.column = "Gene_ID",
#'   sc.filt.genes.cluster = FALSE
#' )
#' SDDLS <- genMixedCellProp(
#'   object = SDDLS,
#'   cell.ID.column = "Cell_ID",
#'   cell.type.column = "Cell_Type",
#'   num.sim.spots = 50,
#'   train.freq.cells = 2/3,
#'   train.freq.spots = 2/3,
#'   verbose = TRUE
#' )
#' SDDLS <- simMixedProfiles(SDDLS)
#' SDDLS <- trainDeconvModel(
#'   object = SDDLS,
#'   batch.size = 12,
#'   num.epochs = 5
#' )
#' }
#'   
trainDeconvModel <- function(
  object,
  type.data.train = "mixed",
  type.data.test = "mixed",
  batch.size = 64,
  num.epochs = 60,
  num.hidden.layers = 2,
  num.units = c(200, 200),
  activation.fun = "relu",
  dropout.rate = 0.25,
  loss = "kullback_leibler_divergence",
  metrics = c("accuracy", "mean_absolute_error",
              "categorical_accuracy"),
  normalize = TRUE,
  scaling = "standardize",
  norm.batch.layers = TRUE,
  custom.model = NULL,
  shuffle = TRUE,
  sc.downsampling = NULL,
  use.generator = FALSE,
  on.the.fly = FALSE,
  agg.function = "AddRawCount",
  threads = 1,
  view.metrics.plot = TRUE,
  verbose = TRUE
) {
  .checkPythonDependencies(alert = "error")
  if (!is(object, "SpatialDDLS")) {
    stop("The provided object is not of SpatialDDLS class")
  } else if (is.null(prob.cell.types(object))) {
    stop("'prob.cell.types' slot is empty")
  } else if (num.epochs <= 2) {
    stop("'num.epochs' argument must be greater than or equal to 2")
  } else if (batch.size < 10) {
    stop("'batch.size' argument must be greater than 10")
  } 
  if (!any(type.data.train %in% c("both", "mixed", "single-cell"))) {
    stop("'type.data.train' argument must be one of the following options: 'both', 'mixed' or 'single-cell'")
  }
  if (!any(type.data.test %in% c("both", "mixed", "single-cell"))) {
    stop("'type.data.test' argument must be one of the following options: 'both', 'mixed' or 'single-cell'")
  }
  # check if data provided is correct regarding on the fly training
  # if (is.null(single.cell.real(object)) && is.null(single.cell.simul(object))) {
  #   stop("At least one single-cell slot must be provided ('single.cell.real' ", 
  #        "or 'single.cell.simul') as trainDeconvModel evaluates ", 
  #        "DNN model on both types of profiles: mixed and single-cell")
  # }
  if (!scaling %in% c("standardize", "rescale", "none")) {
    stop("'scaling' argument must be one of the following options: 'standardize', 'rescale' or 'none'")
  } else {
    if (scaling == "standardize") {
      scaling.fun <- base::scale
    } else if (scaling == "rescale") {
      scaling.fun <- rescale.function
    } else if (scaling == "none") {
      scaling.fun <- function(x) return(x)
    }
  }
  if (!on.the.fly) {
    ## checking if all data are provided
    vec.type.data <- c(type.data.train, type.data.test)
    for (type in seq_along(vec.type.data)) {
      if (verbose & type == 1) message("=== Training and test from stored data")
      text.message <- ifelse(type == 1, "train" , "test")
      if (
        (vec.type.data[type] == "both" && 
         is.null(mixed.profiles(object, type.data = text.message)) ||
         vec.type.data[type] == "both" && 
         (is.null(single.cell.real(object)) && is.null(single.cell.simul(object))))
      ) {
        stop("If `type.data.", text.message, "` = both' is selected, 'mixed.profiles' and at least ",
             "one single cell slot must be provided")
      } else if (vec.type.data[type] == "mixed" && is.null(mixed.profiles(object, type.data = text.message))) {
        stop("If `type.data.", text.message, "` = mixed is selected, 'mixed.profiles' must be provided")
      } 
      
    }  
  } else {
    if (verbose) message("=== Training and test on the fly was selected")
    if ((type.data.train == "both" | type.data.train == "single-cell") && 
        (is.null(single.cell.real(object)) && is.null(single.cell.simul(object)))) {
      stop("If `type.data.train` = both' is selected, at least ",
           "one single cell slot must be provided")
  } 
    ## just in case of on.the.fly = TRUE
    if (!agg.function %in% c("MeanCPM", "AddCPM", "AddRawCount")) {
      stop("'agg.function' must be one of the following options: 'MeanCPM', 'AddCPM', 'AddRawCount'")
    } else {
      if (agg.function == "MeanCPM") {
        .agg.fun <- aggregation.fun.mean.cpm
      } else if (agg.function == "AddCPM") {
        .agg.fun <- aggregation.fun.add.cpm
      } else if (agg.function == "AddRawCount") {
        .agg.fun <- aggregation.fun.add.raw.counts
      }
    }
  }
  if (!is.null(trained.model(object))) {
    warning("'trained.model' slot is not empty. So far, SpatialDDLS",
            " does not support for multiple trained models, so the current model",
            " will be overwritten\n",
            call. = FALSE, immediate. = TRUE)
  }
  # plots in RStudio during training
  if (view.metrics.plot) view.plot <- "auto"
  else view.plot <- 0
  if (verbose) verbose.model <- 1
  else verbose.model <- 0
  prob.matrix.train <- .targetForDNN(
    object = object, combine = type.data.train, 
    shuffle = TRUE, type.data = "train", 
    downsampling = sc.downsampling,
    fly = on.the.fly, verbose = verbose
  )
  prob.matrix.test <- .targetForDNN(
    object = object, combine = type.data.test, 
    shuffle = FALSE, type.data = "test", 
    downsampling = NULL,
    fly = on.the.fly, verbose = verbose
  )
  n.train <- nrow(prob.matrix.train)
  n.test <- nrow(prob.matrix.test)
  # check if the number of samples is compatible with batch.size
  if (n.train < batch.size) {
    stop(
      paste0("The number of samples used for training (", n.train, ") is too ", 
             "small compared with 'batch.size' (", batch.size, "). Please, ", 
             "increase the number of samples or consider reducing 'batch.size'")
    )
  } 
  if (n.test < batch.size) {
    stop(
      paste0("The number of samples used for test (", n.test, ") is too ",
             "small compared with 'batch.size' (", batch.size, "). Please, ",
             "increase the number of samples or consider reducing 'batch.size'")
    )
  }
  if (is.null(custom.model)) {
    if (num.hidden.layers != length(num.units)) {
      stop("The number of hidden layers must be equal to the length of ", 
           "num.units (number of neurons per layer)")
    }
    # check if any argument not provided
    model <- keras_model_sequential(name = "SpatialDDLS")
    for (i in seq(num.hidden.layers)) {
      if (i == 1) {
        model <- model %>% layer_dense(
          units = num.units[i], 
          input_shape = nrow(single.cell.real(object)),
          name = paste0("Dense", i)
        )
      } else {
        model <- model %>% layer_dense(
          units = num.units[i], 
          name = paste0("Dense", i)
        )
      }
      if (norm.batch.layers) {
        model <- model %>% 
          layer_batch_normalization(name = paste0("BatchNormalization", i)) %>%
          layer_activation(activation = activation.fun, 
                           name = paste0("Activation", i)) %>%
          layer_dropout(rate = dropout.rate, name = paste0("Dropout", i))  
      } else {
        model <- model %>% 
          layer_activation(activation = activation.fun, 
                           name = paste0("Activation", i)) %>%
          layer_dropout(rate = dropout.rate, name = paste0("Dropout", i))  
      }
      
    }
    # final layer --> compression and proportions
    if (norm.batch.layers) {
      model <- model %>% layer_dense(
        units = ncol(prob.cell.types(object, "train") %>% prob.matrix()),
        name = paste0("Dense", i + 1)
      ) %>% 
        layer_batch_normalization(name = paste0("BatchNormalization", i + 1)) %>%
        layer_activation(activation = "softmax", name = "ActivationSoftmax")  
    } else {
      model <- model %>% layer_dense(
        units = ncol(prob.cell.types(object, "train") %>% prob.matrix()),
        name = paste0("Dense", i + 1)
      ) %>% layer_activation(activation = "softmax", name = "ActivationSoftmax")  
    }
    
  } else {
    if (!is(custom.model, "keras.engine.sequential.Sequential")) {
      stop("'custom.model' must be a keras.engine.sequential.Sequential object")
    } else if (keras::get_input_shape_at(custom.model$layers[[1]], 1)[[2]] !=
               nrow(single.cell.real(object))) {
      stop("The number of neurons of the first layer must be equal to the ", 
           "number of genes considered by SpatialDDLS object (", 
           nrow(single.cell.real(object))," in this case)")
    } else if (keras::get_output_shape_at(
        custom.model$layers[[length(custom.model$layers)]], 1
      )[[2]] != ncol(prob.cell.types(object, "train") %>% prob.matrix())) {
      stop("The number of neurons of the last layer must be equal to the ", 
           "number of cell types considered by SpatialDDLS object (", 
           ncol(prob.cell.types(object, "train") %>% prob.matrix()), 
           " in this case)")
    } else if (!grepl("'activation': 'softmax'", keras::get_config(custom.model))) {
      stop("In order to get proportions as output, the activation function of the ",
           "last hidden layer must be 'softmax'")
    }
    model <- custom.model
  }
  if (verbose) summary(model)
  # shall I allow to set a custom optimizer?
  model %>% compile(
    loss = loss,
    optimizer = optimizer_adam(),
    metrics = metrics
  )
  if (!is.null(single.cell.simul(object))) {
    suffix.names <- unique(colData(single.cell.simul(object))$suffix)
  } else {
    suffix.names <- "_Simul"
  }
  pattern <- suffix.names
  # set if samples will be generated on the fly
  if (!on.the.fly) {
    .dataForDNN <- .dataForDNN.file
    if (type.data.train == "mixed") {
      checkingClass <- is(
        assay(mixed.profiles(object, type.data = "train")), "HDF5Array"
      )  
    } else if (type.data.train == "single-cell") {
      checkingClass <- is(assay(single.cell.real(object)), "HDF5Array")  
    } else {
      checkingClass <- all(
        is(assay(mixed.profiles(object, type.data = "train")), "HDF5Array"), 
        is(assay(single.cell.real(object)), "HDF5Array")
      )
    }
  } else {
    .dataForDNN <- .dataForDNN.onFly
    checkingClass <- FALSE
  }
  if (verbose) 
    message(paste("\n=== Training DNN with", n.train, "samples:\n"))
  
  if (use.generator | isTRUE(on.the.fly) | checkingClass) {
    gen.train <- .trainGenerator(
      object = object, 
      funGen = .dataForDNN,
      prob.matrix = prob.matrix.train,
      type.data = "train",
      mixing.fun = ifelse(
        any(type.data.train %in% c("mixed", "both")) & on.the.fly == FALSE,
        yes = mixed.profiles(object, "train")@metadata[["mixing.fun"]],
        no = agg.function
      ),
      fun.aggregation = .agg.fun,
      normalize = normalize,
      scaling = scaling.fun,
      batch.size = batch.size,
      combine = type.data.train,
      shuffle = shuffle,
      pattern = pattern,
      min.index = NULL,
      max.index = NULL,
      threads = threads,
      verbose = verbose
    )
    history <- model %>% fit_generator(
      generator = gen.train,
      steps_per_epoch = ceiling(n.train / batch.size),
      epochs = num.epochs,
      verbose = verbose.model,
      view_metrics = view.plot
    )
    if (verbose)
      message(paste0("\n=== Evaluating DNN in test data (", n.test, " samples)"))
    
    # evaluation of the model: set by default, no options?
    gen.test <- .predictGenerator(
      object,
      funGen = .dataForDNN,
      target = TRUE,
      prob.matrix = prob.matrix.test,
      mixing.fun = ifelse(
        any(type.data.test %in% c("mixed", "both")) & on.the.fly == FALSE,
        yes = mixed.profiles(object, "test")@metadata[["mixing.fun"]],
        no = agg.function
      ),
      fun.aggregation = .agg.fun,
      normalize = normalize,
      scaling = scaling.fun,
      batch.size = batch.size,
      pattern = pattern,
      threads = threads,
      verbose = verbose
    )
    test.eval <- model %>% evaluate_generator(
      generator = gen.test,
      steps = ceiling(n.test / batch.size)
    )
    # prediction of test samples
    if (verbose) {
      message(paste0("   - ", names(test.eval), ": ", lapply(test.eval, round, 4),
                     collapse = "\n"))
      message(paste("\n=== Generating prediction results using test data\n"))
    }
    gen.predict <- .predictGenerator(
      object,
      funGen = .dataForDNN,
      target = FALSE,
      prob.matrix = prob.matrix.test,
      mixing.fun = ifelse(
        any(type.data.test %in% c("mixed", "both")) & on.the.fly == FALSE,
        yes = mixed.profiles(object, "test")@metadata[["mixing.fun"]],
        no = agg.function
      ),
      fun.aggregation = .agg.fun,
      normalize = normalize,
      scaling = scaling.fun,
      batch.size = batch.size,
      pattern = pattern,
      threads = threads,
      verbose = verbose
    )
    predict.results <- model %>% predict_generator(
      generator = gen.predict,
      steps = ceiling(n.test / batch.size),
      verbose = verbose.model
    )
    rownames(predict.results) <- rownames(prob.matrix.test)
    colnames(predict.results) <- colnames(prob.matrix.test)  
  } else { # no generators, everything is loaded into memory
    dataTrain <- .dataForDNN.file(
      object = object,
      sel.data = prob.matrix.train,
      pattern = pattern,
      type.data = "train",
      mixing.fun = ifelse(
        any(type.data.train %in% c("mixed", "both")) & on.the.fly == FALSE,
        yes = mixed.profiles(object, "train")@metadata[["mixing.fun"]],
        no = agg.function
      ),
      fun.aggregation = .agg.fun,
      normalize = normalize,
      scaling = scaling.fun,
      threads = threads
    )
    history <- model %>% fit(
      dataTrain,
      prob.matrix.train,
      epochs = num.epochs,
      batch_size = batch.size,
      verbose = verbose.model,
      view_metrics = view.plot
    )
    if (verbose)
      message(paste0("\n=== Evaluating DNN in test data (", n.test, " samples)"))
    
    # evaluation of the model
    dataTest <- .dataForDNN.file(
      object = object,
      sel.data = prob.matrix.test,
      pattern = pattern,
      type.data = "test",
      mixing.fun = ifelse(
        any(type.data.test %in% c("mixed", "both")) & on.the.fly == FALSE,
        yes = mixed.profiles(object, "test")@metadata[["mixing.fun"]],
        no = agg.function
      ),
      fun.aggregation = .agg.fun,
      normalize = normalize,
      scaling = scaling.fun,
      threads = threads
    )
    test.eval <- model %>% evaluate(
      dataTest,
      prob.matrix.test
    )
    # prediction of test samples
    if (verbose) {
      message(paste0("   - ", names(test.eval), ": ", lapply(test.eval, round, 4),
                     collapse = "\n"))
      message(paste("\n=== Generating prediction results using test data\n"))
    }
    predict.results <- model %>% predict(
      dataTest,
      verbose = verbose.model
    )
    rownames(predict.results) <- rownames(prob.matrix.test)
    colnames(predict.results) <- colnames(prob.matrix.test)  
  }
  ## creating DeconvDLModel
  network.object <- new(
    Class = "DeconvDLModel",
    model = model,
    training.history = history,
    test.metrics = test.eval,
    test.pred = predict.results,
    cell.types = colnames(prob.matrix.test),
    features = rownames(single.cell.real(object))
  )
  trained.model(object) <- network.object
  if (verbose) message("DONE")
  
  return(object)
}

.trainGenerator <- function(
  object,
  funGen,
  prob.matrix,
  type.data,
  mixing.fun,
  fun.aggregation,
  normalize,
  scaling,
  batch.size,
  combine,
  shuffle,
  pattern,
  min.index,
  max.index,
  threads,
  verbose
) {
  if (!is.null(min.index) && !is.null(max.index)) {
    n.samples <- length(seq(min.index, max.index))
    nb <- min.index
  } else {
    n.samples <- nrow(prob.matrix)
    min.index <- 1
    nb <- 0
  }
  function() {
    data.index <- seq(nb + 1, nb + batch.size)
    nb <<- nb + batch.size
    if (nb > n.samples) {
      data.index <- data.index[data.index <= n.samples]
      fill <- batch.size - length(data.index)
      data.index <- c(data.index, seq(min.index + 1, min.index + fill))
      if (fill <= min.index) nb <<- min.index + 1
      else nb <<- fill
    }
    if (shuffle) {
      shuffling <- sample(seq_along(data.index))
      sel.data <- prob.matrix[data.index, , drop = FALSE]
      counts <- funGen(
        object = object, 
        sel.data = sel.data, 
        pattern = pattern,
        type.data = type.data,
        mixing.fun = mixing.fun,
        fun.aggregation = fun.aggregation,
        normalize = normalize,
        scaling = scaling,
        threads = threads
      )
      return(list( 
        counts[shuffling, ],
        sel.data[shuffling, ]
      ))
    } else {
      sel.data <- prob.matrix[data.index, , drop = FALSE]
      counts <- funGen(
        object = object, 
        sel.data = sel.data, 
        pattern = pattern,
        type.data = type.data,
        mixing.fun = mixing.fun,
        fun.aggregation = fun.aggregation,
        normalize = normalize,
        scaling = scaling,
        threads = threads
      )
      return(list(counts, sel.data))
    }
  }
}

.predictGenerator <- function(
  object,
  funGen,
  prob.matrix,
  target,
  mixing.fun,
  fun.aggregation,
  normalize,
  scaling,
  batch.size,
  pattern,
  threads,
  verbose
) {
  nb <- 0
  n.samples <- nrow(prob.matrix)
  function() {
    data.index <- seq(nb + 1, nb + batch.size)
    nb <<- nb + batch.size
    if (nb > n.samples) {
      data.index <- data.index[data.index <= n.samples]
      nb <<- 0
    }
    sel.data <- prob.matrix[data.index, , drop = FALSE]
    counts <- funGen(
      object = object, 
      sel.data = sel.data, 
      pattern = pattern,
      type.data = "test",
      mixing.fun = mixing.fun,
      fun.aggregation = fun.aggregation,
      normalize = normalize,
      scaling = scaling,
      threads = threads
    )
    if (target) return(list(counts, sel.data))
    else return(list(counts))
  }
}
## in this function, I have to normalize raw counts
.dataForDNN.file <- function(
  object,
  sel.data,
  pattern,
  type.data,
  mixing.fun,
  fun.aggregation,
  normalize,
  scaling,
  threads
) {
  bulk.data <- grepl(pattern = "Spot_", rownames(sel.data))
  if (any(bulk.data)) {
    bulk.samples <-  as.matrix(
      assay(mixed.profiles(object, type.data))[, rownames(sel.data)[bulk.data],
                                           drop = FALSE]
    )
  } 
  if (any(!bulk.data))  {
    if (!is.null(single.cell.simul(object)) && !is.null(single.cell.real(object))) {
      sel.cells <- rownames(sel.data)[!bulk.data]
      sim.cells <- grep(pattern = pattern, sel.cells, value = TRUE)
      real.cells <- grep(pattern = pattern, sel.cells, value = TRUE, invert = TRUE)
      cell.samples.sim <- as.matrix(
        assay(single.cell.simul(object))[, sim.cells, drop = FALSE]
      )
      cell.samples.real <- as.matrix(
        assay(single.cell.real(object))[, real.cells, drop = FALSE]
      )
      cell.samples <- .mergeMatrices(x = cell.samples.real, y = cell.samples.sim) 
    } else if (!is.null(single.cell.real(object)) && 
               is.null(single.cell.simul(object))) {
      sel.cells <- rownames(sel.data)[!bulk.data]
      cell.samples <- as.matrix(
        assay(single.cell.real(object))[, sel.cells, drop = FALSE]
      )
    } else if (!is.null(single.cell.simul(object)) && 
               is.null(single.cell.real(object))) {
      sel.cells <- rownames(sel.data)[!bulk.data]
      cell.samples <- as.matrix(
        assay(single.cell.simul(object))[, sel.cells, drop = FALSE]
      )
    }
  }
  # return final matrix counts
  if (any(bulk.data) && any(!bulk.data)) {
    cell.samples <- log2(.cpmCalculate(x = cell.samples + 1))
    if (normalize & mixing.fun == "AddRawCount") {
      bulk.samples <- log2(.cpmCalculate(x = bulk.samples + 1))
    }
    counts <- cbind(bulk.samples, cell.samples)[, rownames(sel.data), drop = FALSE]
  } else if (any(bulk.data)) {
    if (normalize & mixing.fun == "AddRawCount") {
      bulk.samples <- log2(.cpmCalculate(x = bulk.samples + 1))
    }
    counts <- bulk.samples[, rownames(sel.data), drop = FALSE]
  } else if (any(!bulk.data)) {
    counts <- log2(
      .cpmCalculate(x = cell.samples[, rownames(sel.data), drop = FALSE] + 1)
    )
  }
  return(scaling(t(counts)))
}

.dataForDNN.onFly <- function(
  object,
  sel.data,
  pattern,
  type.data,
  mixing.fun,
  fun.aggregation,
  normalize,
  scaling,
  threads
) {
  bulk.data <- grepl(pattern = "Spot_", rownames(sel.data))
  if (any(bulk.data)) {
    sel.bulk.cells <- prob.cell.types(object, type.data)@cell.names[
      rownames(sel.data)[bulk.data], , drop = FALSE]
    bulk.samples <- apply(
      X = sel.bulk.cells,
      MARGIN = 1,
      FUN = .setBulk,
      object = object,
      pattern = pattern,
      fun.pseudobulk = fun.aggregation
    )
  } 
  if (any(!bulk.data))  {
    if (!is.null(single.cell.simul(object)) && !is.null(single.cell.real(object))) {
      sel.cells <- rownames(sel.data)[!bulk.data]
      sim.cells <- grep(pattern = pattern, sel.cells, value = TRUE)
      real.cells <- grep(pattern = pattern, sel.cells, value = TRUE, invert = TRUE)
      cell.samples.sim <- as.matrix(
        assay(single.cell.simul(object))[, sim.cells, drop = FALSE]
      )
      cell.samples.real <- as.matrix(
        assay(single.cell.real(object))[, real.cells, drop = FALSE]
      )
      cell.samples <- .mergeMatrices(x = cell.samples.real, y = cell.samples.sim) 
    } else if (!is.null(single.cell.real(object)) && 
               is.null(single.cell.simul(object))) {
      sel.cells <- rownames(sel.data)[!bulk.data]
      cell.samples <- as.matrix(
        assay(single.cell.real(object))[, sel.cells, drop = FALSE]
      )
    } else if (!is.null(single.cell.simul(object)) && 
               is.null(single.cell.real(object))) {
      sel.cells <- rownames(sel.data)[!bulk.data]
      cell.samples <- as.matrix(
        assay(single.cell.simul(object))[, sel.cells, drop = FALSE]
      )
    }
  }
  # return final matrix counts
  if (any(bulk.data) && any(!bulk.data)) {
    cell.samples <- log2(.cpmCalculate(x = cell.samples + 1))
    if (normalize & mixing.fun == "AddRawCount") {
      bulk.samples <- log2(.cpmCalculate(x = bulk.samples + 1))
    }
    counts <- cbind(bulk.samples, cell.samples)[, rownames(sel.data), drop = FALSE]
  } else if (any(bulk.data)) {
    if (normalize & mixing.fun == "AddRawCount") {
      bulk.samples <- log2(.cpmCalculate(x = bulk.samples + 1))
    }
    counts <- bulk.samples[, rownames(sel.data), drop = FALSE]
  } else if (any(!bulk.data)) {
    counts <- log2(
      .cpmCalculate(x = cell.samples[, rownames(sel.data), drop = FALSE] + 1)
    )
  }
  return(scaling(t(counts)))
}

.downsampleCells <- function(
  object, 
  downsampling, 
  type.data, 
  verbose  
) {
  listcells <- prob.cell.types(object, type.data) %>% set.list()
  if (all(downsampling > sapply(listcells, length))) {
    warning(
      "Downsampling is not possible: all cell types have less than ", 
      downsampling, " cells", call. = FALSE, immediate. = TRUE
    )
  } else {
    vecdown <- sapply(listcells, length) > downsampling
    if (any(vecdown)) {
      if (verbose)
        message(
          "=== Downsampling to ", downsampling, " cells for: ", 
          paste(names(vecdown)[vecdown], collapse = ", ")
        )
      for (i in names(vecdown)[vecdown]) {
        listcells[[i]] <- sample(x = listcells[[i]], size = downsampling)  
      }
    } else {
      if (verbose)
        message("=== No cell types with # of cells greater than ", downsampling)
    }  
  }
  return(listcells)
}

.targetForDNN <- function(
  object, 
  combine,
  type.data,
  downsampling,
  fly,
  shuffle,
  verbose
) {
  if (combine == "both") {
    if (!is.null(downsampling)) {
      listcells <- .downsampleCells(
        object = object, downsampling = downsampling, 
        type.data = type.data, verbose = verbose
      )    
    } else {
      listcells <- prob.cell.types(object, type.data) %>% set.list()
    }
    tpsm <- matrix(
      unlist(sapply(
        X = names(listcells),
        FUN = function (x, l) {
          v <- rep(0, length(l))
          names(v) <- names(l)
          v[x] <- 1
          return(rep(v, length(l[[x]])))
        }, 
        l = listcells
      )), 
      ncol = length(listcells), 
      byrow = TRUE,
      dimnames = list(unlist(listcells), names(listcells))
    )
    allCellTypes <- colnames(prob.cell.types(object, type.data) %>% prob.matrix())
    if (!all(allCellTypes %in% colnames(tpsm))) {
      lackTypes <- allCellTypes[!allCellTypes %in% colnames(tpsm)]
      lackMatrix <- matrix(
        0, ncol = length(lackTypes), nrow = nrow(tpsm), 
        dimnames = list(rownames(tpsm), lackTypes)
      )
      tpsm <- cbind(tpsm, lackMatrix)
    }
    tpsm <- tpsm[, colnames(prob.cell.types(object, type.data) %>% 
                              prob.matrix())]

    if (fly) {
      probs.matrix <- rbind(
        tpsm, prob.cell.types(object, type.data) %>% prob.matrix() / 100
      )  
      rownames(probs.matrix) <- c(
        rownames(tpsm), 
        rownames(prob.cell.types(object, type.data) %>% prob.matrix())
      )
    } else {
      tpsm <- tpsm[sample(nrow(tpsm)), ]
      probs.matrix <- prob.cell.types(object, type.data)@prob.matrix[
        colnames(mixed.profiles(object, type.data)), ] / 100
      if (nrow(tpsm) > nrow(probs.matrix)) {
        probs.matrix <- .mergePropsSort(m.small = probs.matrix, m.big = tpsm)
      } else if (nrow(tpsm) <= nrow(probs.matrix)) {
        probs.matrix <- .mergePropsSort(m.small = tpsm, m.big = probs.matrix)
      }
    }
  } else if (combine == "mixed") {
    if (verbose) message("    Using only simulated mixed samples\n")
    if (fly) {
      probs.matrix <- prob.cell.types(object, type.data) %>% prob.matrix() / 100
    } else {
      probs.matrix <- prob.cell.types(object, type.data)@prob.matrix[
        colnames(mixed.profiles(object, type.data)), ] / 100
    }
  } else if (combine == "single-cell") {
    if (verbose) message("    Using only single-cell samples\n")
    if (!is.null(downsampling)) {
      listcells <- .downsampleCells(
        object = object, downsampling = downsampling, 
        type.data = type.data, verbose = verbose
      )    
    } else {
      listcells <- prob.cell.types(object, type.data) %>% set.list()
    }
    
    probs.matrix <- matrix(
      unlist(sapply(
        X = names(listcells),
        FUN = function (x, l) {
          v <- rep(0, length(l))
          names(v) <- names(l)
          v[x] <- 1
          return(rep(v, length(l[[x]])))
        }, l = listcells
      )), ncol = length(listcells), 
      byrow = TRUE,
      dimnames = list(unlist(listcells), names(listcells))
    )
    allCellTypes <- colnames(prob.cell.types(object, type.data) %>% prob.matrix())
    if (!any(allCellTypes %in% colnames(probs.matrix))) {
      lackTypes <- allCellTypes[!allCellTypes %in% colnames(probs.matrix)]
      lackMatrix <- matrix(
        0, ncol = length(lackTypes), nrow = nrow(probs.matrix), 
        dimnames = list(rownames(probs.matrix), lackTypes)
      )
      probs.matrix <- cbind(probs.matrix, lackMatrix)
    }
    probs.matrix <- probs.matrix[
      , colnames(prob.cell.types(object, type.data) %>% prob.matrix())]
  }
  # shuffle only if train on the fly
  if (shuffle && !fly) {
    return(probs.matrix[sample(nrow(probs.matrix)), ])
  } else {
    return(probs.matrix)
  }
}

.mergeMatrices <- function(x, y) {
  genes.out <- setdiff(rownames(x), rownames(y))
  if (identical(genes.out, character(0))) {
    return(cbind(x, y))
  } else {
    zero.genes <- matrix(0, nrow = length(genes.out), ncol = ncol(y), 
                         dimnames = list(genes.out, NULL))
    return(cbind(x, rbind(y, zero.genes)[rownames(x), , drop = FALSE]))  
  }
}

.mergePropsSort <- function(m.small, m.big) {
  nrow.small <- nrow(m.small) 
  nrow.big <- nrow(m.big)
  rows.sel.small <- sort(sample(x = nrow.small + nrow.big, size = nrow.small))
  rows.sel.big <- setdiff(seq(nrow.small + nrow.big), rows.sel.small)
  samples.names <- character()
  samples.names[rows.sel.small] <- rownames(m.small)
  samples.names[rows.sel.big] <- rownames(m.big)
  m.new <- matrix(0L, nrow = nrow.small + nrow.big, ncol = ncol(m.big),
                  dimnames = list(samples.names, colnames(m.big)))
  m.new[rows.sel.small, ] <- m.small
  m.new[rows.sel.big, ] <- m.big
  return(m.new)
}

.simplifySet <- function(vec, index, set) {
  summ <- sum(vec[index])
  # vec <- vec[-c(index)]
  names.vec <- names(vec)
  vec <- c(vec, summ)
  names(vec) <- c(names.vec, set)
  return(vec)
}

.simplifyMajority <- function(vec, index) {
  maxim <- which.max(vec[index])
  summ <- sum(vec[index])
  vec[index[-maxim]] <- 0
  vec[index[maxim]] <- summ
  return(vec)
}

.simplifySetGeneral <- function(results, simplify.set) {
  cell.types <- colnames(results)
  if (any(unlist(lapply(X = simplify.set, FUN = function(x) length(x) < 2))))
    stop("The minimum number of cell types for simplifying is two")
  if (is.null(names(simplify.set)) ||
      (length(names(simplify.set)) != length(simplify.set))) {
    stop("Each element in the list must contain the corresponding new class as name")
  } 
  # else if (length(unique(names(simplify.set))) == length(simplify.set)) {
  #   stop("There are not duplicated names to aggregate results. Items of the list ", 
  #        "must have duplicated names under which to aggregate the results")
  # }
  # check that elements are correct
  lapply(
    X = simplify.set, 
    FUN = function(x, types) {
      if (!all(x %in% types)) {
        stop("Some elements in 'simplify.set' are not present among the cell types ",
             "considered by the model")
      } 
    }, 
    types = cell.types
  )
  index <- lapply(
    X = simplify.set, FUN = function(x) unique(which(colnames(results) %in% x))
  )
  if (any(duplicated(unlist(index)))) {
    stop("'simplify.set' presents duplicated cell types. Please, provide ", 
         "only unique cell types among those considered by the model")
  }
  # for more than 1 subset
  indexNamesUnique <- unique(names(index))
  for (n in indexNamesUnique) {
    results <- t(apply(
      X = results,
      FUN = .simplifySet,
      MARGIN = 1,
      index = index[[n]],
      set = n
    ))
  }
  results <- results[, -unlist(index)]
  return(results)
}

.simplifyMajorityGeneral <- function(results, simplify.majority) {
  cell.types <- colnames(results)
  # check if cell.types provided are correct
  if (any(unlist(lapply(X = simplify.majority, FUN = function(x) length(x) < 2))))
    stop("The minimum number of cell types for simplifying is two")
  lapply(
    X = simplify.majority, 
    FUN = function(x, types) {
      if (!all(x %in% types)) {
        stop("Some elements in 'simplify.majority' are not present between the cell types ",
             "considered by the model")
      } 
    }, 
    types = cell.types
  )
  index <- lapply(
    X = simplify.majority, 
    FUN = function(x) unique(which(colnames(results) %in% x))
  ) 
  if (any(duplicated(unlist(index)))) {
    stop("'simplify.majority' has duplicated cell types. Please, provide ", 
         "only unique cell types among those considered by the model")
  }
  for (n in index) {
    results <- t(apply(
      X = results,
      FUN = .simplifyMajority,
      MARGIN = 1,
      index = n
    ))  
  }
  return(results)
}


################################################################################
############### Deconvolution of spatial transcriptomics data ##################
################################################################################

#' Deconvolute spatial transcriptomics data using trained model
#'
#' Deconvolute spatial transcriptomics data using the trained model in
#' the \code{\linkS4class{SpatialDDLS}} object. The trained model is used 
#' to predict cell proportions of two mirrored transcriptional profiles: 
#' \itemize{ \item 'Intrinsic' profiles: transcriptional profiles of each spot 
#' in the ST dataset. \item 'Extrinsic' profiles: profiles simulated from the 
#' surrounding spots of each spot.} After prediction, cell proportions
#' from the intrinsic profiles (intrinsic cell proportions) are regularized 
#' based on the similarity between intrinsic and extrinsic profiles in order
#' to maintain spatial consistency. This approach leverages both transcriptional 
#' and spatial information. For more details, see Mañanes et al., 2023 and the
#' Details section. 
#' 
#' The deconvolution process involves two main steps: predicting cell 
#' proportions based on transcriptome using the trained neural network model, 
#' and regularization of cell proportions based on the spatial location of each 
#' spot. In the regularization step, a mirrored version of each spot is 
#' simulated based on its N-nearest spots. We refer to these profiles as 
#' 'extrinsic' profiles, whereas the transcriptional profiles of each spot are 
#' called 'intrinsic' profiles. Extrinsic profiles are used to regularize 
#' predictions based on intrinsic profiles. The rationale is that spots 
#' surrounded by transcriptionally similar spots should have similar cell 
#' compositions, and therefore predicted proportions can be smoothed to preserve
#' their spatial consistency. On the other hand, spots surrounded by dissimilar 
#' spots cannot be predicted by their neighbors, and thus they can only be 
#' predicted by their own transcriptional profiles likely due to presenting very 
#' specific cell compositions. 
#' 
#' Regarding the working os \pkg{SpatialDDLS}: first, extrinsic profiles are 
#' simulated based on the N-nearest spots for each spot by summing their 
#' transcriptomes. Distances between extrinsic and intrinsic profiles of each 
#' spot are calculated so that similar/dissimilar spots are identified. These 
#' two sets of transcriptional profiles are used as input for the trained neural
#' network model, and according to the calculated distances, a weighted mean 
#' between the predicted proportions for each spot is calculated. Spots with 
#' distances between intrinsic and extrinsic profiles greater than 
#' \code{alpha.cutoff} are not regularized, whereas spots with distances less 
#' than \code{alpha.cutoff} contribute to the weighted mean. Weights are 
#' calculated by rescaling distances less than \code{alpha.cutoff} between 0 
#' and 0.5, so that the maximum extent to which a extrinsic profile can 
#' modified the predictions based on intrinsic profiles is 0.5 (a regular 
#' mean). For more details, see Mañanes et al., 2023. 
#'
#' This function requires a \code{\linkS4class{SpatialDDLS}} object with a
#' trained deep neural network model (\code{\link{trained.model}} slot, and the
#' spatial transcriptomics datasets to be deconvoluted in the
#' \code{spatial.experiments} slot. See \code{?\link{createSpatialDDLSobject}}
#' or \code{?\link{loadSTProfiles}} for more details.
#'
#' @param object \code{\linkS4class{SpatialDDLS}} object with
#'   \code{trained.model} and \code{spatial.experiments} slots.
#' @param index.st Name or index of the dataset/slide stored in the
#'   \code{SpatialDDLS} object (\code{spatial.experiments} slot) to be
#'   deconvolute. If missing, all datasets will be deconvoluted.
#' @param normalize Normalize data (logCPM) before deconvolution (\code{TRUE} by
#'   default). 
#' @param scaling How to scale data before training. Options include
#'   \code{"standardize"} (values are centered around the mean with a unit
#'   standard deviation) or \code{"rescale"} (values are shifted and rescaled so
#'   that they end up ranging between 0 and 1). If \code{normalize = FALSE},
#'   data are not scaled.
#' @param k.spots Number of nearest spots considered for each spot during 
#'   regularization and simulation of extrinsic transcriptional profiles. The 
#'   greater, the smoother the regularization will be (4 by default). 
#' @param pca.space Whether to use PCA space to calculate distances between 
#'   intrinsic and extrinsic transcriptional profiles (\code{TRUE} by default). 
#' @param fast.pca Whether using the \pkg{irlba} implementation. If \code{TRUE},
#'   the number of PCs used is defined by the \code{} parameter. If 
#'   \code{FALSE}, the PCA implementation from the \pkg{stats} R package is 
#'   used instead (\code{TRUE} by default).
#' @param pcs.num Number of PCs used to calculate distances if 
#'   \code{fast.pca == TRUE} (50 by default). 
#' @param pca.var Threshold of explained 
#'   variance (between 0.2 and 1) used to choose the number of PCs used if 
#'   \code{pca.space == TRUE} and \code{fast.pca == FALSE} (0.8 by default). 
#' @param metric Metric used to measure distance/similarity between intrinsic 
#'   and extrinsic transcriptional profiles. It may be \code{'euclidean'}, 
#'   \code{'cosine'} or \code{'pearson'} (\code{'euclidean'} by default).
#' @param alpha.cutoff Minimum distance for regularization. 
#'   It may be \code{'mean'} (spots with transcriptional distances shorter than 
#'   the mean distance of the dataset will be modified) or \code{'quantile'} 
#'   (spots with transcriptional distances shorter than the 
#'   \code{alpha.quantile} quantile are used). \code{'mean'} by default.
#' @param alpha.quantile Quantile used if \code{alpha.cutoff == 'quantile'}. 
#'   0.5 by default.  
#' @param simplify.set List specifying which cell types should be compressed
#'   into a new label with the name of the list item. See examples for details.
#'   If provided, results are stored in a list with \code{'raw'} and
#'   \code{'simpli.set'} elements.
#' @param simplify.majority List specifying which cell types should be
#'   compressed into the cell type with the highest proportion in each spot.
#'   Unlike \code{simplify.set}, no new labels are created. If provided, results
#'   are stored in a list with \code{'raw'} and \code{'simpli.majority'}
#'   elements.
#' @param use.generator Boolean indicating whether to use generators for
#'   prediction (\code{FALSE} by default).
#' @param batch.size Number of samples per batch. Only when \code{use.generator
#'   = TRUE}.
#' @param verbose Show informative messages during the execution.
#'
#' @return \code{\linkS4class{SpatialDDLS}} object with a \code{deconv.spots}
#'   slot. The output is a list containing 'Regularized', 'Intrinsic' and 
#'   'Extrinsic' deconvoluted cell proportions, 'Distances' between intrinsic 
#'   and extrinsic transcriptional profiles, and 'Weight.factors' with the 
#'   final weights used to regularize intrinsic cell proportions. If 
#'   \code{simplify.set} and/or \code{simplify.majority} are provided, 
#'   the \code{deconv.spots} slot will contain a list with raw and simplified 
#'   results.
#'
#' @export
#'
#' @seealso \code{\link{trainDeconvModel}} \code{\linkS4class{SpatialDDLS}}
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(
#'     counts = matrix(
#'      rpois(30, lambda = 5), nrow = 15, ncol = 20,
#'       dimnames = list(paste0("Gene", seq(15)), paste0("RHC", seq(20)))
#'     )
#'   ),
#'   colData = data.frame(
#'     Cell_ID = paste0("RHC", seq(20)),
#'     Cell_Type = sample(x = paste0("CellType", seq(6)), size = 20,
#'                        replace = TRUE)
#'   ),
#'   rowData = data.frame(
#'     Gene_ID = paste0("Gene", seq(15))
#'   )
#' )
#' SDDLS <- createSpatialDDLSobject(
#'   sc.data = sce,
#'   sc.cell.ID.column = "Cell_ID",
#'   sc.gene.ID.column = "Gene_ID",
#'   sc.filt.genes.cluster = FALSE
#' )
#' SDDLS <- genMixedCellProp(
#'   object = SDDLS,
#'   cell.ID.column = "Cell_ID",
#'   cell.type.column = "Cell_Type",
#'   num.sim.spots = 50,
#'   train.freq.cells = 2/3,
#'   train.freq.spots = 2/3,
#'   verbose = TRUE
#' ) 
#' SDDLS <- simMixedProfiles(SDDLS)
#' # training of SDDLS model
#' SDDLS <- trainDeconvModel(
#'   object = SDDLS,
#'   batch.size = 15,
#'   num.epochs = 5
#' )
#' # simulating spatial data
#' ngenes <- sample(3:40, size = 1)
#' ncells <- sample(10:40, size = 1)
#' counts <- matrix(
#'   rpois(ngenes * ncells, lambda = 5), ncol = ncells,
#'   dimnames = list(paste0("Gene", seq(ngenes)), paste0("Spot", seq(ncells)))
#' )
#' coordinates <- matrix(
#'   rep(c(1, 2), ncells), ncol = 2
#' )
#' st <- SpatialExperiment::SpatialExperiment(
#'   assays = list(counts = as.matrix(counts)),
#'   rowData = data.frame(Gene_ID = paste0("Gene", seq(ngenes))),
#'   colData = data.frame(Cell_ID = paste0("Spot", seq(ncells))),
#'   spatialCoords = coordinates
#' )
#' SDDLS <- loadSTProfiles(
#'   object = SDDLS,
#'   st.data = st,
#'   st.spot.ID.column = "Cell_ID",
#'   st.gene.ID.column = "Gene_ID"
#' )
#' # simplify arguments
#' simplify <- list(CellGroup1 = c("CellType1", "CellType2", "CellType4"),
#'                  CellGroup2 = c("CellType3", "CellType5"))
#' SDDLS <- deconvSpatialDDLS(
#'   object = SDDLS,
#'   index.st = 1,
#'   simplify.set = simplify, 
#'   simplify.majority = simplify
#' )
#' }
#' @references Mañanes, D., Rivero-García, I., Jimenez-Carretero, D., 
#'   Torres, M., Sancho, D., Torroja, C., Sánchez-Cabo, F. (2023). SpatialDDLS: 
#'   An R package to deconvolute spatial transcriptomics data using neural 
#'   networks. biorxiv. doi: \doi{10.1101/2023.08.31.555677}.
#'   
deconvSpatialDDLS <- function(
  object,
  index.st,
  normalize = TRUE,
  scaling = "standardize",
  k.spots = 4, 
  pca.space = TRUE,
  fast.pca = TRUE,
  pcs.num = 50,
  pca.var = 0.8,
  metric = "euclidean",
  alpha.cutoff = "mean",
  alpha.quantile = 0.5,
  simplify.set = NULL,
  simplify.majority = NULL,
  use.generator = FALSE,
  batch.size = 64,
  verbose = TRUE
) {
  # check if python dependencies are covered
  .checkPythonDependencies(alert = "error")
  if (!is(object, "SpatialDDLS")) {
    stop("The provided object is not of class SpatialDDLS")
  } else if (is.null(trained.model(object))) {
    stop("There is no trained model in SpatialDDLS object")
  } else if (is.null(spatial.experiments(object))) {
    stop("No spatial transcriptomics data provided. See ?createSpatialDDLSobject or ?loadSTProfiles")
  }
  if (!scaling %in% c("standardize", "rescale", "none")) {
    stop("'scaling' argument must be one of the following options: 'standardize' or 'rescale'")
  } else {
    if (scaling == "standardize") {
      scaling.fun <- base::scale
    } else if (scaling == "rescale") {
      scaling.fun <- rescale.function
    } else if (scaling == "none") {
      scaling.fun <- function(x) return(x)
    }
  }
  ## checking parameters for regularization
  if (pca.var > 1 | pca.var <= 0.2) {
    stop("pca.var must be a doublet between 0.2 and 1")
  }
  if (!any(metric %in% c("euclidean", "cosine", "pearson"))) {
    stop("metric only can be one of the following options: 'euclidean', 'cosine', 'pearson'")
  }
  if (!any(alpha.cutoff %in% c("mean", "quantile"))) {
    stop("alpha.cutoff only can be 'mean' or 'quantile'")
  }
    
  # checking if DNN model is json format or compiled
  if (is.list(trained.model(object)@model)) {
    model.comp <- .loadModelFromJSON(trained.model(object))
    trained.model(object) <- model.comp
  }
  if (missing(index.st)) {
    if (verbose) {
      message(
        "   No 'index.st' provided. Deconvoluting all SpatialExperiment objects ", 
        "contained in the `spatial.experiments` slot\n"
      ) 
    }
    index.st <- seq_along(spatial.experiments(object))
  } else {
    if (is.character(index.st) & !is.null(names(spatial.experiments(object)))) {
      ## check if all index.st are present in the slot
      stopifnot(
        "`index.st` contains elements not present in spatial.experiments slot " = index.st %in% 
          names(spatial.experiments(object))
      )
    }
  }
  if (is.character(index.st)) {
    namesList <- index.st
  } else {
    namesList <- names(spatial.experiments(object))[index.st]
  }
  deconv.st <- lapply(
    X = index.st,
    FUN = function(index) {
      deconv.counts <- assay(spatial.experiments(object, index))
      results <-.deconvCore(
        deconv.counts = deconv.counts,
        model = trained.model(object),
        batch.size = batch.size,
        normalize = normalize,
        scaling = scaling.fun,
        use.generator = use.generator,
        verbose = verbose
      )
      ## spatial consistency: add here the function (with a conditional I guess)
      list.res.post <- .spatialRegularization(
        spatial.experiments(object, index.st = index), 
        prop.intrinsic = results,
        model = trained.model(object),
        batch.size = batch.size,
        normalize = normalize,
        scaling = scaling.fun,
        use.generator = use.generator,
        k = k.spots, 
        pca = pca.space,
        fast.pca = fast.pca,
        pcs.num = pcs.num,
        pca.var = pca.var,
        metric = metric,
        alpha.cutoff = alpha.cutoff,
        alpha.quantile = alpha.quantile,
        verbose = verbose
      )
      list.res.post[["Instrinsic.Proportions"]] <- results
      ## simplfy: not very used to be honest, let's keep it just for convenience 
      if (!is.null(simplify.set) || !is.null(simplify.majority)) {
        if (verbose) message("=== Note that only regularized proportions will be simplified")
        
        list.res <- list(raw = list.res.post)
        if (!is.null(simplify.set)) {
          if (!is(simplify.set, "list")) {
            stop("'simplify.set' must be a list in which each element is a ", 
                 "cell type considered by the model")
          }
          results.set <- .simplifySetGeneral(
            results = list.res.post[["Regularized"]],
            simplify.set = simplify.set
          )
          list.res[["simpli.set"]] <- results.set
        }
        if (!is.null(simplify.majority)) {
          if (!is(simplify.majority, "list")) {
            stop("'simplify.majority' must be a list in which each element is a ", 
                 "cell type considered by the model")
          }
          results.maj <- .simplifyMajorityGeneral(
            results = list.res.post[["Regularized"]],
            simplify.majority = simplify.majority
          )
          list.res[["simpli.majority"]] <- results.maj
        }
        return(list.res)
      } 
      return(list.res.post)
    }
  ) %>% setNames(namesList)
  deconv.spots(object) <- deconv.st
  if (verbose) message("DONE")
  
  return(object)
}

.deconvCore <- function(
  deconv.counts,
  model,
  batch.size,
  normalize,
  scaling,
  use.generator,
  verbose
) {
  if (is.null(rownames(deconv.counts))) {
    stop("The given matrix does not have column names. You must provide a matrix",
         " with feature names in the same notation used in training data")
  }
  # this can do it more elegantly and efficiently
  # filtering features missing in training data
  filter.features <- rownames(deconv.counts) %in% features(model)
  deconv.counts <- deconv.counts[filter.features, ]
  fill.features <- !features(model) %in% rownames(deconv.counts)
  if (any(fill.features)) {
    m.new <- matrix(0L, nrow = sum(fill.features), ncol = ncol(deconv.counts))
    rownames(m.new) <- features(model)[fill.features]
    deconv.counts <- rbind(deconv.counts, m.new)
    deconv.counts <- deconv.counts[features(model), ]  
  }
  if (verbose) {
    if (sum(!filter.features) != 0) {
      message(
        paste(
          "=== Filtering out", sum(!filter.features),
          "features in data that are not present in trained model\n"
        )
      )  
    }
    if (sum(fill.features) != 0) {
      message(
        paste(
          "=== Setting", sum(fill.features), 
          "features that are not present in trained model to zero\n"
        )
      ) 
    }
  }
  if (normalize) {
    if (verbose) message("=== Normalizing data (LogCPM)\n")
    
    deconv.counts <- log2(.cpmCalculate(x = deconv.counts + 1))
    deconv.counts <- scaling(t(deconv.counts))
  } else {
    deconv.counts <- as.matrix(deconv.counts)
    deconv.counts <- scaling(t(deconv.counts))
  }
  
  if (verbose) {
    verbose.model <- 1
    message("=== Predicting cell type proportions\n")
  } else {
    verbose.model <- 0
  }
  dnn.model <- model(model)
  if (use.generator) {
    deconv.generator <- .predictDeconvDataGenerator(
      data = deconv.counts,
      model = model,
      batch.size = batch.size
    )
    results <- dnn.model %>% predict_generator(
      generator = deconv.generator,
      steps = ceiling(nrow(deconv.counts) / batch.size),
      verbose = verbose.model
    )  
  } else {
    results <- dnn.model %>% predict(
      deconv.counts,
      verbose = verbose.model
    )  
  }
  if (!is.null(rownames(deconv.counts))) {
    rownames.deconv <- rownames(deconv.counts)
  } else {
    rownames.deconv <- seq(1, nrow(deconv.counts))
  }
  rownames(results) <- rownames.deconv
  colnames(results) <- cell.types(model)

  return(results)
}

.predictDeconvDataGenerator <- function(
  data,
  model,
  batch.size
) {
  nb <- 0
  n.samples <- nrow(data)
  n.features <- length(features(model))
  n.classes <- length(cell.types(model))
  function() {
    data.index <- seq(nb + 1, nb + batch.size)
    nb <<- nb + batch.size
    if (nb > n.samples) {
      data.index <- data.index[data.index <= n.samples]
      nb <<- 0
    }
    return(list(matrix(data[data.index, ],
                       ncol = n.features,
                       nrow = length(data.index))))
  }
}

.spatialRegularization <- function(
    st.data, 
    prop.intrinsic,
    model,
    scaling,
    batch.size = batch.size,
    normalize = normalize,
    use.generator = use.generator,
    k = k, 
    pca = pca,
    fast.pca = fast.pca,
    pcs.num = pcs.num,
    pca.var = pca.var,
    metric = metric,
    alpha.cutoff = alpha.cutoff,
    alpha.quantile = alpha.quantile,
    verbose = verbose
) {
  ## calculating kNN spots per spot
  distances <- as.matrix(dist(SpatialExperiment::spatialCoords(st.data)))
  knn.indexes <- FNN::get.knn(distances, k = k)[[1]]
  # print(knn.indexes)
  rownames(knn.indexes) <- rownames(SpatialExperiment::spatialCoords(st.data))
  
  ## generating mixed transcriptional profiles using kNN spots
  summ.nn.expr <- apply(
    X = knn.indexes, MARGIN = 1, 
    FUN = \(idx) Matrix::rowSums(SummarizedExperiment::assay(st.data)[, idx])
  )
  colnames(summ.nn.expr) <- colnames(SummarizedExperiment::assay(st.data))
  ## prediction of cell proportions on extrinsic signals
  prop.nn <-.deconvCore(
    deconv.counts = summ.nn.expr,
    model = model,
    batch.size = batch.size,
    normalize = normalize,
    scaling = scaling,
    use.generator = use.generator,
    verbose = FALSE
  )
  ## normalizing both data in order to calculate distances 
  logcpm.nn.expr <- log2(.cpmCalculate(summ.nn.expr + 1))
  logcpm.spot.expr <- log2(
    .cpmCalculate(
      SummarizedExperiment::assay(st.data) + 1
    )
  )
  ## calculating distances between actual and contextual matrices
  if (isTRUE(pca)) {
    if (verbose) message("\n=== Calculating distances in PCA space")
    
    logcpm.nn.expr.2 <- logcpm.nn.expr
    colnames(logcpm.nn.expr.2) <- paste0(colnames(logcpm.nn.expr), "_NN")
    
    if (fast.pca) {
      if (!requireNamespace("irlba", quietly = TRUE)) {
        warning(
          paste(
            "irlba is required but not available. Running", 
            "PCA implementation from the stats R package"
          ), call. = FALSE, immediate. = TRUE
        )
      }
      if (verbose) message("\n=== Calculating ", pcs.num, " PCs\n")
      
      log.cpm.final <- cbind(logcpm.spot.expr, logcpm.nn.expr.2)
      if (pcs.num >= min(ncol(log.cpm.final), nrow(log.cpm.final))) {
        pcs.num <- ceiling(min(nrow(log.cpm.final), ncol(log.cpm.final)) * 0.8)
      }
      pca.space <- irlba::irlba(A = t(x = log.cpm.final), nv = pcs.num)
      pca.nn <- pca.space[[2]][grepl("_NN$", colnames(log.cpm.final)), ]
      pca.ori <- pca.space[[2]][!grepl("_NN$", colnames(log.cpm.final)), ]
    } else {
      pca.space <- stats::prcomp(
        t(cbind(logcpm.spot.expr, logcpm.nn.expr.2)), center = TRUE, scale. = TRUE
      )
      eigs <- pca.space$sdev^2
      n.pcs <- sum(cumsum(eigs / sum(eigs)) <= pca.var)
      if (verbose) 
        message(paste0("        - Using ", n.pcs, " PCs (variance cutoff: ", pca.var, ")\n"))
        
      pca.nn <- pca.space$x[grep("_NN$", rownames(pca.space$x)), 1:n.pcs]
      pca.ori <- pca.space$x[grep("_NN$", rownames(pca.space$x), invert = T), 1:n.pcs]
    }
    ## distances
    dist.mm <- .distCalc(x = pca.ori, y = pca.nn, metric = metric)
  } else {
    if (verbose) message("\n=== Calculating distances in transcriptome space\n")
    
    dist.mm <- .distCalc(x = logcpm.spot.expr, y = logcpm.nn.expr, metric = metric)
  }
  ## calculating alpha regularization based on distances
  if (verbose) message("=== Calculating alpha factors based on distances\n")
  
  if (alpha.cutoff == "mean") {
    dist.mod <- as.vector(dist.mm) <= mean(as.vector(dist.mm))
    dist.rescaled <- rescale.function.2(
      as.vector(dist.mm)[dist.mod], 
      to = c(0, 0.5)
    )  
  } else if (alpha.cutoff == "quantile") {
    if (alpha.quantile > 0.99 | alpha.quantile < 0.01) {
      stop("alpha.quantile only can be a number between 0.99 and 0.01")
    }
    dist.mod <- as.vector(dist.mm) <= stats::quantile(
      x = as.vector(dist.mm), probs = alpha.quantile
    )
    dist.rescaled <- rescale.function.2(
      as.vector(dist.mm)[dist.mod], 
      to = c(0, 0.5)
    ) 
  }
  alpha.extrin <- ifelse(test = dist.mod, yes = dist.rescaled, no = 0)
  alpha.intric <- 1 - alpha.extrin 
  
  ## probably, this part can be solved using matrix multiplication, check in the future
  regulaized.prop <- matrix(0, nrow = nrow(prop.intrinsic), ncol = ncol(prop.intrinsic))
  for (spot in seq(nrow(prop.intrinsic))) {
    regulaized.prop[spot, ] <- sapply(
      X = colnames(prop.intrinsic), 
      FUN = \(x) {
        stats::weighted.mean(
          c(prop.intrinsic[spot, x], prop.nn[spot, x]), 
          w = c(alpha.intric[spot], alpha.extrin[spot])
        )
      }
    )  
  }
  colnames(regulaized.prop) <- colnames(prop.intrinsic)
  rownames(regulaized.prop) <- rownames(prop.intrinsic)
  return(
    list(
      regulaized.prop, prop.intrinsic, prop.nn, as.vector(dist.mm), alpha.extrin
    ) %>% setNames(
        c("Regularized", "Intrinsic", 
          "Extrinsic", "Distances", "Alpha.Factors")
      )
  )
}

.distCalc <- function(x, y, metric) {
  dist.mm <- matrix(0, ncol = 1, nrow = nrow(x))
  if (metric == "euclidean") {
    for (i in seq(nrow(x))) {
      dist.mm[i, 1] <- stats::dist(
        t(cbind(x[i, ], y[i, ])), method = "euclidean"
      )
    }  
  } else if (metric == "cosine") {
    for (i in seq(nrow(x))) {
      if (!requireNamespace("lsa", quietly = TRUE)) {
        stop("lsa R package is required but not available")
      }
      dist.mm[i, 1] <- lsa::cosine(x = x[i, ], y = y[i, ])
    } 
    dist.mm <- 1 - dist.mm ## from similarity to distance
  } else if (metric == "pearson") {
    for (i in seq(nrow(x))) {
      dist.mm[i, 1] <- stats::cor(
        x = x[i, ], y = y[i, ], method = "pearson"
      )
    } 
    dist.mm[dist.mm < 0] <- 0
    dist.mm <- 1 - dist.mm 
  }
  
  return(dist.mm)
}

rescale.function.2 <- function(
    x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE)
) {
  (x - from[1]) / diff(from) * diff(to) + to[1]
}

rescale.function <- function(x) {
  apply(X = x, MARGIN = 2, FUN = function(x) (x - min(x)) / (max(x) - min(x)))
}

.loadModelFromJSON <- function(object) {
  model.list <- model(object)
  model.comp <- model_from_json(model.list[[1]])
  model.comp <- set_weights(model.comp, model.list[[2]])
  model(object) <- model.comp
  return(object)
}

.saveModelToJSON <- function(object) {
  model.comp <- model(object)
  model.json <- model_to_json(model.comp)
  weights <- get_weights(model.comp)
  model(object) <- list(model.json, weights)
  return(object)
}
