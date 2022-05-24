# AdaptGauss2D
Adapting two-dimensional Gaussian Mixture Models on data with the help of an interactive visualization and a manual controllable expectation maximization algorithm.

### R shiny app
The AdaptGauss2D app is a GUI made for an intuitive approach for two-dimensional gaussian mixture modeling (GMM). There is also an approach for the one-dimensional case, called AdaptGauss (https://cran.r-project.org/web/packages/AdaptGauss/AdaptGauss.pdf). The output of the app can be used as input. Be aware, that the principal component axis and the angle are not required, since this information is within the covariance matrices and is only provided in the output as a quick access of the results. The essence of AdaptGauss2D lies in the coarse manual fitting as a first step and a fine automatized fitting with an EM algorithm as the second step. That way the user can communicate a rough expectation of the model which is then optimized. In other words, the user provides a rough desired initialization and the automaton of AdaptGauss2D finds a local optimum for that. AdaptGauss2D yields an initialized solution for a GMM as proposal and start for the user.

### Practical impact
Experience showed that Expectations Maximization (EM) algorithms do not necessarily yield the "best" options. In the one- or two-dimensional data case EM algorithms oftentimes create models which seem crooked and off when visualized. Therefore, we created an interactive visualization tool, with which you can adapt models to the data both with automatic optimization and with user intervention. The visualization tool and a performance measure (here: root mean square deviation RMSD) enable a feedback for the user. The rich tools offered by the GUI support the user in an intuitive way to cope with the model fitting task.

### Big Data
The initial computation on the app start require some time for large datasets, but the efficient sampling techniques used for visualization support a comparably fast manual control during app use.

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
