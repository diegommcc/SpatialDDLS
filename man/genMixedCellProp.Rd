% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simMixedSpots.R
\name{genMixedCellProp}
\alias{genMixedCellProp}
\title{Generate training and test cell type composition matrices}
\usage{
genMixedCellProp(
  object,
  cell.ID.column,
  cell.type.column,
  num.sim.spots,
  n.cells = 50,
  train.freq.cells = 3/4,
  train.freq.spots = 3/4,
  proportion.method = c(0, 0, 1),
  prob.sparity = 1,
  min.zero.prop = NULL,
  balanced.type.cells = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{object}{\code{\linkS4class{SpatialDDLS}} object with
\code{single.cell.real} slot and, optionally, with \code{single.cell.simul}
slot.}

\item{cell.ID.column}{Name or column number corresponding to cell names in
cells metadata.}

\item{cell.type.column}{Name or column number corresponding to cell types in
cells metadata.}

\item{num.sim.spots}{Number of mixed profiles to be simulated. It is
recommended to adjust this number according to the number of available
single-cell profiles.}

\item{n.cells}{Specifies the number of cells to be randomly selected and
combined to generate the simulated mixed profiles. By default, it is set to
50 It controls the level of noise present in the simulated data, as it
determines how many single-cell profiles will be combined to produce each
spot.}

\item{train.freq.cells}{Proportion of cells used to simulate training mixed
transcriptional profiles (3/4 by default).}

\item{train.freq.spots}{Proportion of mixed transcriptional profiles to be
used for training, relative to the total number of simulated spots
(\code{num.sim.spots}). The default value is 3/4.}

\item{proportion.method}{Vector with three elements that controls the
proportion of simulated proportions generated by each method: random
sampling of a Dirichlet distribution, "pure" spots (1 cell type), and spots
generated from a random sampling of a Dirichlet distribution but with a
specified number of different cell types (determined by
\code{min.zero.prop}), respectively. By default, all samples are generated
by the last method.}

\item{prob.sparity}{It only affects the proportions generated by the first
method (Dirichlet distribution). It determines the probability of having
missing cell types in each simulated spot, as opposed to a mixture of all
cell types. A higher value for this parameter will result in more sparse
simulated samples.}

\item{min.zero.prop}{This parameter controls the minimum number of cell types
that will be absent in each simulated spot. If \code{NULL} (by default),
this value will be half of the total number of different cell types, but
increasing it will result in more spots composed of fewer cell types. This
helps to create more sparse proportions and cover a wider range of
situations during model training.}

\item{balanced.type.cells}{Boolean indicating whether training and test cells
will be split in a balanced way considering cell types (\code{TRUE} by
default).}

\item{verbose}{Show informative messages during the execution (\code{TRUE} by
default).}
}
\value{
A \code{\linkS4class{SpatialDDLS}} object with \code{prob.cell.types}
slot containing a \code{list} with two \code{\linkS4class{PropCellTypes}}
objects (training and test). For more information about the structure of
this class, see \code{?\linkS4class{PropCellTypes}}.
}
\description{
Generate training and test cell type composition matrices for the simulation
of mixed transcriptional profiles with known cell composition using
single-cell expression profiles. The resulting
\code{\linkS4class{PropCellTypes}} object will contain all the information
needed to simulate new mixed transcriptional profiles. Note this function
does not simulate the mixed profiles, this task is performed by the
\code{\link{simMixedProfiles}} or \code{\link{trainDeconvModel}} functions
(see Documentation).
}
\details{
First, the single-cell profiles are randomly divided into two subsets, with
2/3 of the data for training and 1/3 for testing. The default setting for
this ratio can be changed using the \code{train.freq.cells} parameter. Next,
a total of \code{num.sim.spots} mixed proportions are simulated using a
Dirichlet distribution. This simulation takes into account the probability of
missing cell types in each spot, which can be adjusted using the
\code{prob.sparity} parameter. For each mixed sample, \code{n.cells}
single-cell profiles are randomly selected and combined to generate the
simulated mixed sample. In addition to the Dirichlet-based proportions, pure
spots (containing only one cell type) and spots containing a specified number
of different cell types (determined by the \code{min.zero.prop} parameter)
are also generated in order to cover situations with only a few cell types
present. The proportion of simulated spots generated by each method can be
controlled using the \code{proportion.method} parameter. To visualize the
distribution of cell type proportions generated by each method, the
\code{\link{showProbPlot}} function can be used.
}
\examples{
set.seed(123)
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(
    counts = matrix(
      rpois(100, lambda = 5), nrow = 40, ncol = 30,
      dimnames = list(paste0("Gene", seq(40)), paste0("RHC", seq(30)))
    )
  ),
  colData = data.frame(
    Cell_ID = paste0("RHC", seq(30)),
    Cell_Type = sample(x = paste0("CellType", seq(4)), size = 30,
                       replace = TRUE)
  ),
  rowData = data.frame(
    Gene_ID = paste0("Gene", seq(40))
  )
)

SDDLS <- createSpatialDDLSobject(
  sc.data = sce,
  sc.cell.ID.column = "Cell_ID",
  sc.gene.ID.column = "Gene_ID",
  sc.filt.genes.cluster = FALSE,
  project = "Simul_example"
)
SDDLS <- genMixedCellProp(
  object = SDDLS,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_Type",
  num.sim.spots = 10,
  train.freq.cells = 2/3,
  train.freq.spots = 2/3,
  verbose = TRUE
)
  
}
\seealso{
\code{\link{simMixedProfiles}} \code{\linkS4class{PropCellTypes}}
}
