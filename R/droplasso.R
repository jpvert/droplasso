#' Fit a droplasso model
#'
#' Fit a dropout lasso (droplasso) model.
#'
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is an
#'   observation vector.
#' @param y Response variable.
#' @param family Response type. \code{family="gaussian"} (default) for least squares regression, \code{family="binomial"} for logistic regression
#' @param keep_prob The probability that each element is kept (default: \code{0.5})
#' @param nlambda The number of \code{lambda} values (default: \code{100})
#' @param lambda.min.ratio Smallest value for \code{lambda}, as a fraction of
#'   \code{lambda.max}, the (data derived) entry value (i.e. the smallest value
#'   for which all coefficients are zero). The default depends on the sample
#'   size \code{nobs} relative to the number of variables \code{nvars}. If
#'   \code{nobs > nvars}, the default is \code{0.0001}, close to zero.  If
#'   \code{nobs < nvars}, the default is \code{0.01}.
#' @param lambda The sequence of regularization parameters. By default,
#'   \code{lambda=NULL} lets the function estimate a good sequence by itself.
#' @param init Initial model to start optimization (default: zero vector).
#' @param gamma0 Initial value of the learning rate (default: \code{1})
#' @param decay Learning rate decay (default: \code{0.01})
#' @param n_passes Number of passes over each example of the data on average (default: \code{1000})
#' @param minibatch_size Batch size (default: \code{nobs})
#' @return A vector of coefficients of size \code{nvars} that solves the droplasso problem. 
#' @details The optimization problem is solved with a stochastic proximal gradient descent algorithm, using mini-batches of size \code{minibatch_size}, and a learning rate decaying as \code{gamma0/(1+decay*t)}, where \code{t} is the number of mini-batches processed.
#' @examples
#' #create data:
#' nobs = 100
#' nvars = 5
#' x = matrix(rnorm(nobs*nvars),nrow=nobs)
#' b = c(1,1,0,0,0)
#' p = 1/(1+exp(-x%*%b))
#' y = p>0.5
#' # Fit a lasso model (no dropout)
#' droplasso(x, y, family="binomial", lambda=0.1, keep_prob=1)
#' # Fit a dropout model (no lasso)
#' droplasso(x, y, family="binomial", lambda=0, keep_prob=0.5)
#' # Fit a dropout lasso model
#' droplasso(x, y, family="binomial", lambda=0.1, keep_prob=0.5)
#' @export
droplasso <- function(x, y, family=c("gaussian","binomial"), keep_prob=0.5, nlambda=10, lambda.min.ratio = ifelse(nobs<nvars,0.01,0.0001), lambda=NULL, init=matrix(0,nrow=ncol(x)), gamma0=1, decay=0.01, n_passes=1000, minibatch_size=nrow(x))
{
  this.fcall=match.call()
  family = match.arg(family)
  if (keep_prob > 1) {
    warning("keep_prob >1; set to 1")
    keep_prob = 1
  }
  if (keep_prob < 0) {
    warning("keep_prob<0; set to 0")
    keep_prob = 0
  }
  keep_prob = as.double(keep_prob)
  
  if (is.null(lambda)) {
    lambda_max <- lambda_max(x,y,family)
    lambda_min=lambda_max*lambda.min.ratio
    lambda=lambda_max*(lambda_min/lambda_max)^(c(0:(nlambda-1))/(nlambda-1))
    skipfirstlambda <- TRUE
  } else {
    nlambda=length(lambda)
    skipfirstlambda <- FALSE
  }
  
  if (skipfirstlambda) {
    beta <- matrix(0,nrow=ncol(x))
  } else {
    beta <- droplassoC(x, y, family, keep_prob, lambda[1], init, gamma0, decay, n_passes, minibatch_size)
  }
  
  if (nlambda>1) {
    for (ilambda in seq(2,nlambda)) {
      beta <- cbind(beta, droplassoC(x, y, family, keep_prob, lambda[ilambda], beta[,ilambda-1], gamma0, decay, n_passes, minibatch_size))
    }
  }
  
  res=list(beta=beta, lambda=lambda, nzero=apply(beta!=0, 2, sum), call=this.fcall)
  class(res) <- 'droplasso'
  return(res)
}


#' Compute the smallest lambda for which 0 is solution (i.e, the largest lambda
#' of the regularization path). Check the paper for the derivation. Note that 
#' in the case of logistic loss, we assume the class coded as 0/1; and in the
#' case of gaussian loss, we consider 1/2*RSS
lambda_max <- function(x, y, family=c("gaussian","binomial")) {
  family = match.arg(family)
  offset <- switch(family,
                   gaussian=0,
                   binomial=0.5)
  n <- dim(x)[1]
  s <- max(abs(t(x) %*% (y-offset)))/n
  return(s)
}
