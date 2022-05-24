# AdaptGauss2D
Adapting two-dimensional Gaussian Mixture Models on data with the help of an interactive visualization and a manual controllable expectation maximization algorithm.

# Practical impact
Experience showed that Expectations Maximization (EM) algorithms do not necessarily yield the "best" options. In the one- or two-dimensional data case EM algorithms oftentimes create models which seem crooked and off when visualized. Therefore, we created an interactive visualization tool, with which you can adapt models to the data both with automatic optimization and with user intervention. The visualization tool and a performance measure (here: root mean square deviation RMSD) enable a feedback for the user. The rich tools offered by the GUI support the user in an intuitive way to cope with the model fitting task.


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
