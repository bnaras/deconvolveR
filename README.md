# Empirical Bayes Deconvolution

An unknown prior density $g(\theta)$ has yielded (unobservable) $\Theta_1, \Theta_2,\ldots,\Theta_N$, and each $\Theta_i$ produces
an observation $X_i$ from an exponential family. `deconvolveR` is an R package for estimating prior distribution $g(\theta)$ from the data
using Empirical Bayes deconvolution.

The current package is still under construction but will soon appear
on CRAN along with a manuscript. Meanwhile, you can reproduce many
examples by installing the package in R thus:

```
devtools::install_github("bnaras/deconvolveR")
```


