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


# SpatialDDLS 1.0.2 (2024-03-21)

* Included two parameters to control python and tensorflow versions in the installTFpython function.
* The installTFpython function now installs python 3.8 by default. 

# SpatialDDLS 1.0.3 (2024-10-28)

* SingleCellExperiment adds a saveRDS S3 generics in the new BioC version. To avoid incompatibilities, it is not completely loaded anymore. 
