#' @import methods
NULL

#' SpatialDDLS: an R package to deconvolute spatial transcriptomics data using
#' deep neural networks
#'
#' \pkg{SpatialDDLS} is an R package that provides a neural network-based
#' solution for cell type deconvolution of spatial transcriptomics data. The
#' package takes advantage of single-cell RNA sequencing (scRNA-seq) data to
#' simulate mixed transcriptional profiles with known cell composition and train
#' fully-connected neural networks to predict the cell type composition of
#' spatial transcriptomics spots. The resulting trained models can be applied to
#' new spatial transcriptomics data to predict cell type proportions, allowing
#' for more accurate cell type identification and characterization of
#' spatially-resolved transcriptomic data. Finally, predictions are forced to keep 
#' spatial consistency through a process we refer to as spatial regularization. 
#' Overall, \pkg{SpatialDDLS} is a powerful tool for cell type deconvolution in 
#' spatial transcriptomics data, providing a reliable, fast and flexible 
#' solution for researchers in the field. See Mañanes et al. (2023) 
#' (\doi{10.1101/2023.08.31.555677}) and some examples 
#' (\url{https://diegommcc.github.io/SpatialDDLS/}) for more details.
#'
#' @docType package
#' @name SpatialDDLS-Rpackage
NULL
#> NULL
