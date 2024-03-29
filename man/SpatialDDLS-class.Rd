% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{SpatialDDLS-class}
\alias{SpatialDDLS-class}
\alias{SpatialDDLS}
\title{The SpatialDDLS Class}
\description{
The \code{\linkS4class{SpatialDDLS}} object is the core of the
\pkg{SpatialDDLS} package. This object stores different intermediate data
needed for the construction of new deconvolution models, the spatial
transcriptomics profiles to be deconvoluted, and the predicted cell type
proportions.
}
\details{
This object uses other classes to store different types of data generated
during the workflow: \itemize{ \item
\code{\linkS4class{SingleCellExperiment}} class for single-cell RNA-Seq data
storage, using sparse matrix from the \pkg{Matrix} package
(\code{\linkS4class{dgCMatrix}} class) or \code{HDF5Array} class in case of
using HDF5 files as back-end (see below for more information). \item
\code{\linkS4class{SpatialExperiment}} class for spatial transcriptomics data
storage. \item \code{\linkS4class{ZinbModel}} class with estimated parameters
for the simulation of new single-cell profiles. \item
\code{\linkS4class{SummarizedExperiment}} class for simulated mixed
transcriptional profiles storage. \item \code{\linkS4class{PropCellTypes}}
class for composition cell type matrices. See
\code{?\linkS4class{PropCellTypes}} for details. \item
\code{\linkS4class{DeconvDLModel}} class to store information related to
deep neural network models. See \code{?\linkS4class{DeconvDLModel}} for
details. }

In order to provide a way to work with large amounts of data in
RAM-constrained machines, we provide the possibility of using HDF5 files as
back-end to store count matrices of both real and simulated single-cell
profiles by using the \pkg{HDF5Array} and \pkg{DelayedArray} classes from the
homonymous packages.
}
\section{Slots}{

\describe{
\item{\code{single.cell.real}}{Real single-cell data stored in a
\code{SingleCellExperiment} object. The count matrix is stored either as
\code{\linkS4class{dgCMatrix}} or \code{HDF5Array} objects.}

\item{\code{spatial.experiments}}{List of \code{\linkS4class{SpatialExperiment}}
objects to be deconvoluted.}

\item{\code{zinb.params}}{\code{\linkS4class{ZinbModel}} object with estimated
parameters for the simulation of new single-cell expression profiles.}

\item{\code{single.cell.simul}}{Simulated single-cell expression profiles using the
ZINB-WaVE model.}

\item{\code{prob.cell.types}}{\code{\linkS4class{PropCellTypes}} class with cell
composition matrices built for the simulation of mixed transcriptional
profiles with known cell composition.}

\item{\code{mixed.profiles}}{List of simulated train and test mixed transcriptional
profiles. Each entry is a \code{\linkS4class{SummarizedExperiment}} object.
Count matrices can be stored as \code{HDF5Array} objects using HDF5 files
as back-end in case of RAM limitations.}

\item{\code{trained.model}}{\code{\linkS4class{DeconvDLModel}} object with
information related to the deconvolution model. See
\code{?\linkS4class{DeconvDLModel}} for more details.}

\item{\code{deconv.spots}}{Deconvolution results. It consists of a list where each
element corresponds to the results for each
\code{\linkS4class{SpatialExperiment}} object contained in the
\code{spatial.experiments} slot.}

\item{\code{project}}{Name of the project.}

\item{\code{version}}{Version of \pkg{SpatialDDLS} this object was built under.}
}}

