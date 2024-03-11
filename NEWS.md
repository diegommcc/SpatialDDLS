# SpatialDDLS 0.1.0 (2023-04-08)

* Added a `NEWS.md` file to track changes.


# SpatialDDLS 0.2.0 (2023-10-04)

* Added a vignette explaining HDF5 file usage (hdf5Backend.Rmd vignette). 
* Mixed transcriptional profiles are now stored as raw counts rather than 
normalized values in order to make calculations more transparent. 
* Scale factor for normalization can be chosen (10e3 is now the default option). 


# SpatialDDLS 1.0.0 (2023-12-05)

* Regularization of predicted cell proportions incorporated. Functions and classes relying on the deconvSpatialDDLS function have been modified. 
* Added a set of functions for clustering analysis based on predicted cell proportions (spatialClustering.R file).
* Added a module for neural network interpretation based on the vanilla gradient algorithm (interGradientsDL.R file).
* Changes in default parameters and vignette updated.


# SpatialDDLS 1.0.1 (2024-02-07)

* Change in HDF5 file usage: new version of the HDF5Array package does not support "for.use" argument.
* In createSpatialDDLSobject, included the sc.log.FC parameter to optionally choose if filtering genes according to logFC. 
