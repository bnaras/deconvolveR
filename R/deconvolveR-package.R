#' R package for Empirical Bayes \eqn{g}-modeling using exponential families.
#'
#' \code{deconvolveR} is a package for Empirical Bayes Deconvolution
#' and Estimation. A friendly introduction is provided in the JSS
#' paper reference below and this package includes a vignette
#' containing a number of examples.
#'
#' @docType package
#' @name deconvolveR-package
#'
#' @references Bradley Efron. Empirical Bayes Deconvolution
#'     Estimates. Biometrika 103(1), 1-20, ISSN
#'     0006-3444. doi:10.1093/biomet/asv068.
#'     \url{http://biomet.oxfordjournals.org/content/103/1/1.full.pdf+html}
#' @references Bradley Efron and Trevor Hastie. Computer Age
#'     Statistical Inference.  Cambridge University Press. ISBN
#'     978-1-1-7-14989-2. Chapter 21.
#' @references Balasubramanian Narasimhan and Bradley
#'     Efron. deconvolveR: A G-Modeling Program for Deconvolution and
#'     Empirical Bayes Estimation. doi:10.18637/jss.v094.i11
#'
NULL

#' Shakespeare word counts in the entire canon: 14,376 distinct words
#' appeared exactly once, 4343 words appeared twice etc.
#'
#' @references Bradley Efron and Ronald Thisted. Estimating the number of
#' unseen species: How many words did Shakespeare know? Biometrika, Vol 63(3),
#' doi:10.1093/biomet/63.3.435.
#'
#' @name bardWordCount
#' @usage data(bardWordCount)
#' @docType data
#' @keywords data
#'
NULL

#' Intestinal surgery data involving 844 cancer patients. The data consists of
#' pairs (\eqn{n_i}, \eqn{s_i}) where \eqn{n_i} is the number of satellites removed
#' and \eqn{s_i} is the number of satellites found to be malignant.
#'
#' @references Gholami, et. al. Number of Lymph Nodes Removed and Survival after
#' Gastric Cancer Resection: An Analysis from the US Gastric Cancer Collaborative.
#' J Am Coll Surg. 2015 Aug;221(2):291-9. doi: 10.1016/j.jamcollsurg.2015.04.024.
#'
#' @name surg
#' @docType data
#' @usage data(surg)
#' @keywords data
#'
NULL

#' A set of \eqn{\Theta} values that have a bimodal distribution for testing
#'
#' @name disjointTheta
#' @docType data
#' @usage data(disjointTheta)
#' @keywords data
#'
NULL

