# deconvolveR

## Empirical Bayes Deconvolution

An unknown prior density \\g(\theta)\\ has yielded (unobservable)
\\\Theta_1, \Theta_2,\ldots,\Theta_N\\, and each \\\Theta_i\\ produces
an observation \\X_i\\ from an exponential family. `deconvolveR` is an R
package for estimating prior distribution \\g(\theta)\\ from the data
using Empirical Bayes deconvolution.

Details and examples may be found in the paper by [Narasimhan and Efron,
2020](https://dx.doi.org/10.18637/jss.v094.i11). A vignette with further
examples is also provided.
