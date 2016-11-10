#' A function to compute Empirical Bayes estimates using deconvolution
#'
#' @importFrom splines ns
#' @importFrom stats nlm dpois dnorm dbinom
#' @param tau a vector of (implicitly m) discrete support points for \eqn{\theta}
#' @param X the vector of sample values: a vector of counts for Poisson, a vector of z-scores for Normal,
#'        a 2-d matrix with rows consisting of pairs, trial size, \eqn{n_i}, and number of successes,
#'        \eqn{X_i}, for Binomial
#' @param Q the Q matrix, implies y and P are supplied as well
#' @param P the P matrix, implies Q and y are supplied as well
#' @param n the number of bins in the discretization (default 40)
#' @param c0 the regularization parameter (default 1)
#' @param family the exponential family, one of \code{c("Poisson", "Normal", "Binomial")}
#'        with \code{"Poisson"}, the default
#' @param ignoreZero if the zero values should be ignored (default = \code{TRUE}). Applies to
#'        Poisson only.
#' @param deltaAt the theta value where a delta function is desired (default \code{NULL}). Applies to
#'        Normal only if non-null
#' @param scale if the Q matrix should be scaled so that the spline basis has mean 0 and columns
#'        sum of squares to be one, (default \code{TRUE})
#' @param pDegree the degree of the splines to use (default 5)
#' @param aStart the starting values for the non-linear optimization, default is a vector of 1s
#' @param ... further args to function \code{nlm}
#' @return a list of 9 items consisting of
#' \item{mle}{the maximum likelihood estimate \eqn{\hat{\alpha}}}
#' \item{Q}{the m by p matrix Q}
#' \item{P}{the n by m matrix P}
#' \item{R}{the ratio of artificial to genuine information as per the reference below}
#' \item{cov}{the covariance matrix for the mle}
#' \item{cov.g}{the covariance matrix for the \eqn{g}}
#' \item{mat}{an m by 6 matrix containing columns for \eqn{theta}, \eqn{g}, std. error of \eqn{g},
#'            \eqn{G} (the cdf of {g}), std. error of \eqn{G}, and the bias of \eqn{g}}
#' \item{loglik}{the negative log-likelihood function for the data taking a \eqn{p}-vector argument}
#' \item{statsFunction}{a function to compute the statistics returned above}
#' @export
#' @references Bradley Efron. Empirical Bayes Deconvolution Estimates. Biometrika 103(1), 1-20,
#' ISSN 0006-3444. doi:10.1093/biomet/asv068.
#' \url{http://biomet.oxfordjournals.org/content/103/1/1.full.pdf+html}
#' @references Bradley Efron and Trevor Hastie. Computer Age Statistical Inference.
#' Cambridge University Press. ISBN 978-1-1-7-14989-2. Chapter 21.
#' @examples
#' set.seed(238923) ## for reproducibility
#' N <- 1000
#' theta <- rchisq(N,  df = 10)
#' nSIM <- 50
#' data <- sapply(seq_len(nSIM), function(x) rpois(n = N, lambda = theta))
#' tau <- seq(1, 32)
#' results <- apply(data, 2,
#'                  function(x) deconv(theta = tau, x = x, ignoreZero = FALSE,
#' stats <- sapply(results, function(x) x$stats[, "g"])
#' mean <- apply(stats, 1, mean)
#' sd <- apply(stats, 1, sd)
#' data.frame(theta = tau, gTheta = 100 * dchisq(tau, df = 10),
#'            mean = 100 * mean, sd = 100 * sd, cv = sd / mean)
#'            c0 = 1))
deconv <- function(tau, X, y, Q, P, n = 40, family = c("Poisson", "Normal", "Binomial"),
                   ignoreZero = TRUE, deltaAt = NULL, c0 = 1,
                   scale = TRUE,
                   pDegree = 5,
                   aStart = 1.0, ...) {

    if (missing(Q) && missing(P)) {
        family <- match.arg(family)

        m <- length(tau)

        if (family == "Poisson") {
            if (ignoreZero) {
                x0 <- seq_len(n)
                ## Accounting for truncation at 0
                P <- sapply(tau, function(x) dpois(x = x0, lambda = x) / (1 - exp(-x)) )
            } else {
                x0 <- seq.int(from = 0, to = n - 1)
                P <- sapply(tau, function(x) dpois(x = x0, lambda = x) )
            }

            if (missing(y)) {
                y <- sapply(x0, function(i) sum(X == i))
            }

            ##Q <- cbind(sqrt((m - 1) / (m)), scale(splines::ns(tau, pDegree)))
            Q <- cbind(1.0, scale(splines::ns(tau, pDegree), center = TRUE, scale = FALSE))
            Q <- apply(Q, 2, function(w) w / sqrt(sum(w * w)))

        } else if (family == "Normal") {
            ## x is the z scores
            r <- round(range(X), digits = 1)
            xBin <- seq(from = r[1], to = r[2], length.out = n)
            xBinDropFirst <- xBin[-1]
            xBinDropLast <- xBin[-length(xBin)]

            P <- sapply(tau,
                        function(x) pnorm(q = xBinDropFirst, mean = x) -
                                    pnorm(q = xBinDropLast, mean = x)
                        )
            intervals <- findInterval(X, vec = xBin)
            y <- sapply(seq_len(n-1), function(w) sum(intervals == w))

            if (scale) {
                Q1 <- scale(splines::ns(tau, pDegree), center = TRUE, scale = FALSE)
                Q1 <- apply(Q1, 2, function(w) w / sqrt(sum(w * w)))
            }
            if (!is.null(deltaAt)) {
                I0 <- as.numeric(abs(tau - deltaAt) < 1e-10)
                Q <- cbind(I0, Q1)
            }  else {
                Q <- Q1
            }

        } else { ## family == "Binomial"
            Q <- splines::ns(tau, pDegree)
            if (scale) {
                Q <- scale(Q, center = TRUE, scale = FALSE)
                Q <- apply(Q, 2, function(w) w / sqrt(sum(w * w)))
            }

            ## x is 2-d: n_i, s_i
            P <- sapply(tau,
                        function(w)  dbinom(X[, 2], size = X[, 1], prob = w))
            y <- 1

        }
    } else { ## user specified Q or  P
        if (missing(y) || missing(P) || missing(Q)) {
            stop("P, Q, and y must be specified together!")
        }
    }

    p <- ncol(Q)
    pGiven <- length(aStart)
    if (pGiven == 1) {
        aStart <- rep(aStart, p)
    } else {
        if (pGiven != p)
            stop(sprintf("Wrong length (%d) for initial parameter, expecting length (%d)",
                         pGiven, p))
    }

    statsFunction <- function(a) {
        ## Penalized MLE bias and cov statistics
        ## N=1 for "unique" case 6/10/14

        ## See rtest/substitute.progs for my "dot functions"

        g <- as.vector(exp(Q %*% a))
        g <- g / sum(g)
        G <- cumsum(g)
        f <- as.vector(P %*% g)
        yHat <- sum(y) * f
        Pt <- P / f
        W <- g * (t(Pt) - 1 )
        qw <- t(Q) %*% W
        ywq <- (yHat * t(W)) %*% Q
        I1 <- qw %*% ywq   ## Fisher Information Matrix I(\alpha)

        aa <- sqrt(sum(a^2))
        sDot <- c0 * a / aa
        sDotDot <- (c0 / aa) * ( diag(length(a)) - outer(a, a) /aa^2 )

        ## The R value
        R <- sum(diag(sDotDot)) / sum(diag(I1))

        I2 <- solve(I1 + sDotDot)
        bias <- as.vector(-I2 %*% sDot)
        Cov <- I2 %*% (I1 %*% t(I2))
        ## Se <- diag(Cov)^.5
        Dq <- (diag(g) - outer(g, g)) %*% Q
        bias.g <- Dq %*% bias
        Cov.g <- Dq %*% Cov %*% t(Dq)
        se.g <- diag(Cov.g)^.5

        D <- diag(length(tau))
        D[lower.tri(D)] <- 1
        Cov.G <- D %*% (Cov.g %*% t(D))
        se.G <- diag(Cov.G)^.5
        mat <- cbind(tau, g, se.g, G, se.G, bias.g)
        colnames(mat) = c("theta", "g", "SE.g", "G", "SE.G", "Bias.g")
        list(R = R, cov = Cov, cov.g = Cov.g, mat = mat)
    }

    loglik <- function(a) {
        ##-loglikelihood for g-modeling mle calcs 12/20,21/12
        ##  y0=1 for Unique Case

        g <- exp(Q %*% a)
        g <- as.vector(g / sum(g))
        f <- as.vector(P %*% g)
        value <- -sum(y * log(f)) + c0 * sum(a^2)^.5
        Pt <- P / f
        W <- g * (t(Pt) - 1 )
        qw <- t(Q) %*% W
        aa <- sqrt(sum(a^2))
        sDot <- c0 * a / aa
        if (family == "Binomial") {
            attr(value, "gradient") <- -rowSums(qw) + sDot
        } else {
            attr(value, "gradient") <- - (qw %*% y)  + sDot
        }
        ## yHat <- sum(y) * f
        ## ywq <- (yHat * t(W)) %*% Q
        ## I1 <- qw %*% ywq   ## Fisher Information Matrix I(\alpha)
        ## sDotDot <- (c0 / aa) * ( diag(length(a)) - outer(a, a) / aa^2 )
        ## I2 <- solve(I1 + sDotDot)
        ## attr(value, "hessian") <- -I2
        value
    }

    result <- stats::nlm(f = loglik, p = aStart, gradtol = 1e-10, ...)

    mle <- result$estimate
    stats <- statsFunction(mle)
    list(mle = mle,
         Q = Q,
         P = P,
         R = stats$R,
         cov = stats$cov,
         cov.g = stats$cov.g,
         stats = stats$mat,
         loglik = loglik,
         statsFunction = statsFunction)
}






