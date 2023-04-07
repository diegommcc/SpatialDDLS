# SpatialDDLS: development notebook

Here, I write down ideas, future work, or general things I think are worth considering about the package. 

## Ideas

* Saliency maps: as now I'm correctly scaling my data before training, maybe this approach could work. Take one mode and check it out on Python. 
* Smoothing proportions considering spatial features: using gradient descent, I implemented some code able to minimize each proportion independently. This could be interesting if the proper spatial information is taken. 

## Future work

* Both training and test are done by batches and generators: 
    * Batches are needed for training, not for test. Could this be changed?ç
    * The usage of generators make all calculations slower compared with loading everything once. I think I should change this without removing the code, i.e., only use it when a lot of data is provided. 

**About generators** **DONE**

I think that considering how the package is being used, the generators are not needed anymore, just when on.the.fly or HDF5 files are used. Therefore, I should implement alternative ways to train the model. **DONE**

## Things I should compulsory check before release it on CRAN

* Simplify proportions after prediction: (simpli.set and simpli.majority)
    * I think this part is not checked, so probably there are some mistakes that should be taken into account. 
    * Moreover, this part should be check: Idk what happens by default.
* HDF5 files: there must be mistakes because this code has not changed at all. Check it on the test. 
* Change all examples
* Tests
* Change %>% by |>: not possible, it is implemented from R 4.1 
