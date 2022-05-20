# AdaptGauss2D
Adapting two-dimensional Gaussian Mixture Models to data with help of interactive visualization and single expectation maximization steps.

## Installation using Github
Please note, that dependecies have to be installed manually.

```R
remotes::install_github("Mthrun/AdaptGauss2D")
```

## Tutorial Examples

The tutorial with several examples can be found on in the video 
https://www.youtube.com/watch?v=MV7DVEWys_c

If the program was started with “results=AdaptGauss2D(Data)”, the variable results will be a list of the following elements: 
“Means”, “CovarianceMatrices”, “Weights”, “PrincipalComponentAxis”, “Angle” and “Cls”. 
If the given dataset has more than two features, a different pair of features can be chosen during the app runtime and a new analysis can be started.
