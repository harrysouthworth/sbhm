#' @import rjags coda
NULL

runStan <- function(data, nchains, iter, burn, thin) {
    is.stan <- requireNamespace("rstan")
    if (!is.stan) {
        stop("rstan is not available")
    }
    mod <- "data {  int<lower=0> J;
                    real y[J];
                    real<lower=0> sigma[J];
                    real m;
                    real tauScale;
                  }
            parameters {
                    real theta[J];
                    real mu;
                    real<lower=0> tau;
            }
            model {
                    mu ~ normal(m, 10000);          // standard deviation is 10^4
                    tau ~ cauchy(0, tauScale);      // scale
                    theta ~ normal(mu, tau);
                    y ~ normal(theta, sigma);
            }"
    rstan::stan(model_code=mod, data=data,
                chains=nchains, iter=iter, warmup=burn, thin=thin)
}

runJags <- function(data, nchains, iter, burn, thin) {
    ## jags is always available (from the NAMESPACE)
    mod <- "model {
              for(i in 1:J) {
                theta[i] ~ dnorm(mu, tau)
                thes[i] <- pow(sigma[i], -2)
                y[i] ~ dnorm(theta[i], thes[i])
              }
              mu ~ dnorm(m, pow(10000, -2))
              itau ~ dt(0, pow(tauScale, -2), 1) T(0,)
              tau <- pow(itau, -2)
            }"
    sink(tempfile()) # Hide output: suppressMessages won't work
    on.exit(sink())
    res <- jags.model(textConnection(mod),
                      data=data,
                      n.chains=nchains)
    coda.samples(res, n.iter=iter - burn,
                 c("theta", "mu", "tau"),
                 thin=thin)
}

#' Fit a simple hierarchical normal model with MCMC
#' @aliases print.sbhm summary.sbhm plot.sbhm print.summary.sbhm get.sbhm.summary
#' @param y The observed group means.
#' @param s The observed group standard errors.
#' @param m The top level mean. Defaults to \code{m=0}.
#' @param tauScale The scale parameter for the half-Cauchy distribution for the top
#'   level scale. Defaults to twice the range of \code{y}.
#' @param engine Which MCMC engine to use. Defaults to \code{engine="stan"} and
#'   if \code{stan} is not avaiable will try \code{JAGS}. Note that if the engine
#'   is JAGS, a since Markov chain will be run, irrespective of the \code{nchains}
#'   argument.
#' @param nchains The number of chains to run. If not specified and \code{stan} is
#'   available, the function tries to figure out how many cores are available,
#'   then uses all but one of them to run a chain. If JAGS is to be used, it defaults
#'   to running 3 chains in sequence, not parallel.
#' @param iter The number of steps in each Markov chain. Defaults to \code{iter=41000}.
#' @param burn The number of burn-in steps to be discarded from each chain.
#'   Defaults to \code{burn=1000}.
#' @param thin The amount of thinning of each chain to do. Defaults to \code{thin=4}.
#' @param probs The quantiles of the posterior distribution to print. Defaults to
#'   \code{probs=c(.025, .25, .5, .75, .975)}.
#' @param x In \code{plot.sbhm}, an object of class "sbhm".
#' @param xlab In \code{plot.sbhm}, the x-axis label.
#' @param ylab In \code{plot.sbhm}, the y-axis label.
#' @param main In \code{plot.sbhm}, the main title.
#' @param offset In \code{plot.sbhm}, how much to nudge the mle estiamtes upwards and the
#'   sbhm estimates downwards. Defaults to \code{offset=.2}.
#' @param mle.col In \code{plot.sbhm}, the colour of the mles.
#' @param sbhm.col In \code{plot.sbhm}, the colour of the shrunken estimates.
#' @param cex, In \code{plot.sbhm}, character expansion. Defaults to \code{cex=1.5}.
#' @param margins Vector of length 4 giving the margin sizes for \code{plot}. Defaults
#'   to \code{margins=c(5.1, 7.1, 4.1, 2.1)}.
#' @param lwd.step The difference in line widths for the confidence intervals. Defaults
#'   to \code{lwd.step=2}.
#' @param digits In \code{print.sbhm}, the number of digits to round to.
#' @param ... Additional arguments to \code{plot}. Not used.
#' @return An object of class "sbhm" containing the fitted model, \code{probs},
#'   \code{engine}, the data and function call.
#' @details The function requires that you have either rstan or rjags installed.
#'   Because many users of rstan would have no need for rjags, rjags is not
#'   featured as a dependency. Because rstan is not on CRAN and takes a little
#'   effort to install, the user might want to use rjags instead, so rstan is
#'   not featured as a dpeendency.
#' @author Harry Southworth
#' @examples
#' \dontrun{
#' mod <- sbhm(schools$effect, schools$se)
#' plot(mod)
#' # The MERIT-HF study
#' y <- merit$logHR
#' names(y) <- merit$country # used as labels in the plot
#' mod <- sbhm(y=y, s=merit$se, nchains=2)
#' plot(mod)
#' }
#' @export sbhm
sbhm <- function(y, s, m=0, tauScale=NULL, engine=c("stan", "jags"),
                 nchains=NULL, iter=41000, burn=1000, thin=4,
                 probs=c(.025, .25, .5, .75, .975)){
    thecall <- match.call()
    
    engine <- match.arg(engine)

    ## setup the data
    
    J <- length(y)
    if (length(s) != J) {
        stop("y and s have different lengths")
    }
    
    if (missing(tauScale)) {
        tauScale <- 2 * abs(diff(range(y))) # let it fail if there are NAs
    }
    
    d <- list(J=J, y=y, sigma=s, tauScale=tauScale, m=m)
    
    res <-
        switch(engine,
               stan={
                   nchains <- if (is.null(nchains)) {1} else {nchains}
                   runStan(d, nchains, iter, burn, thin)
               },
               jags={
                   nchains <- if (is.null(nchains)) {3} else {nchains}
                   runJags(d, nchains, iter, burn, thin)
               },
               stop("unknown MCMC engine"))
    
    data <- data.frame(y=y, s=s)
    if (!is.null(names(y))) rownames(data) <- names(y)
    
    res <- list(fit=res,
                probs=probs,
                engine=engine,
                data=data,
                call=thecall)

    class(res) <- "sbhm"
    res
}


#' @export get.sbhm.summary
get.sbhm.summary <- function(x){
    switch(x$engine,
           stan={
               res <- rstan::summary(x$fit, probs=x$probs)[[1]]
               ## Drop log-posterior row, various suammry cols
               res[-nrow(res), 4:8]
           },
           jags={
               res <- summary(x$fit, quantiles=x$probs)[[2]]
               ## reorder columns to match stan output
               res <- rbind(res[-c(1:2), ], res[1:2, ])
               ## convert to same parameterization as stan
               res[nrow(res), ] <- rev(sqrt(1/res[nrow(res), ]))
               res
           },
           stop("unknown engine"))
}

#' @describeIn sbhm Print method for sbhm
#' @export
print.sbhm <- function(x, digits=3, ...){
  res <- get.sbhm.summary(x)
  print(res, digits=digits)
  invisible(NULL)
}

#' @export
summary.sbhm <- function(object, ...){
    ## default stan or jags output. Assume the user knows what they're doing
    res <- switch(object$engine,
                  stan=rstan::summary(object$fit, probs=object$probs),
                  jags=summary(object$fit, quantile=object$probs))
    class(res) <- "summary.sbhm"
    res
}

##' @export
print.summary.sbhm <- function(x, ...){
  print(unclass(x))
  invisible(x)
}

#' @describeIn sbhm Plot an sbhm object
#' @export
plot.sbhm <- function(x, xlab="Estimates", ylab="", main="", mle.col="blue", sbhm.col="orange",
                     margins=c(5.1, 7.1, 4.1, 2.1),
                     offset=.2, lwd.step=2, cex=1.5, ...){
  # Resize the margins
  on.exit(oldpar)
  oldpar <- par(no.readonly=TRUE)
  par(mar=margins)
  
  # First get the initial estimates and put intervals on them
  d <- x$data
  n <- length(d$y)
  p <- x$probs
  alpha <- p[p > .5] # actually 1 - alpha/2
  z <- sort(qnorm(alpha))
  # Get confidence intervals
  ci <- t(sapply(1:n, function(X) sort(c(d$y[X] - z * d$s[X], d$y[X], d$y[X] + z * d$s[X]))))
  nci <- (ncol(ci) - 1) / 2 # number of CIs
  lwds <- seq(1, 1 + (nci - 1) * lwd.step, by=lwd.step) # line widths to use in plotting CIs

  ests <- get.sbhm.summary(x)[1:n, ]

  # Estimates from hierarchical model will be shrunk so x-axis scale needs only
  # to contain ci, so plot it first

  yy <- (n:1) - offset
  plot(d$y, yy, xlim=range(ci), pch=16, col=mle.col, xlab=xlab, ylab=ylab, main=main, cex=cex,
       axes=FALSE, ylim=range(c(yy, (n:1) + offset)))
  box()
  axis(1)
  axis(2, at=n:1, labels=rownames(d), las=2)
  for (i in 1:nci)
    segments(x0=ci[, i], x1=ci[, ncol(ci) - i + 1], y0=yy, lwd=lwds[i], col=mle.col)
  yy <- (n:1) + offset
  points(ests[, nci + 1], yy, pch=16, col=sbhm.col, cex=cex)
  for (i in 1:nci)
    segments(x0=ests[, i], x1=ests[, ncol(ests) - i + 1], y0=yy, lwd=lwds[i], col=sbhm.col)

  invisible(NULL)
}
