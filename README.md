# deconvolveR <img src="man/figures/logo.png" align="right" height="139" alt="deconvolveR hex logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/bnaras/deconvolveR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bnaras/deconvolveR/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/deconvolveR)](https://cran.r-project.org/package=deconvolveR)
[![Monthly downloads](https://cranlogs.r-pkg.org/badges/deconvolveR)](https://cran.r-project.org/package=deconvolveR)
[![Total downloads](https://cranlogs.r-pkg.org/badges/grand-total/deconvolveR)](https://cran.r-project.org/package=deconvolveR)
<!-- badges: end -->

## Empirical Bayes Deconvolution

An unknown prior density $g(\theta)$ has yielded (unobservable)
$\Theta_1, \Theta_2,\ldots,\Theta_N$, and each $\Theta_i$ produces an
observation $X_i$ from an exponential family. `deconvolveR` is an R
package for estimating prior distribution $g(\theta)$ from the data
using Empirical Bayes deconvolution.

Details and examples may be found in the paper by [Narasimhan and Efron,
2020](https://dx.doi.org/10.18637/jss.v094.i11). A vignette with further
examples is also provided.
