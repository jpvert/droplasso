#' Fit a droplasso model
#'
#' Fit a dropout lasso (droplasso) model.
#'
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is an
#'   observation vector.
#' @param y Response variable.
#' @param family Response type. \code{family="gaussian"} (default) for least squares regression, \code{family="binomial"} for logistic regression
#' @param keep_prob The probability that each element is kept (default: \code{0.5})
#' @param lambda Regularisation parameter of the l_1 norm (default: \code{0})
#' @param init Initial model to start optimization (default: zero vector).
#' @param gamma0 Initial value of the learning rate (default: \code{1})
#' @param decay Learning rate decay (default: \code{0.01})
#' @param n_passes Number of passes over each example of the data on average (default: \code{1000})
#' @param minibatch_size Batch size (default: \code{nobs})
#' @return A vector of coefficients of size \code{nvars} that solves the droplasso problem. The optimization problem is solved with a stochastic proximal gradient descent algorithm, using mini-batches of size \code{minibatch_size}, and a learning rate decaying as \code{gamma0/(1+decay*t)}, where \code{t} is the number of mini-batches processed.
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
droplasso <- function(x, y, family=c("gaussian","binomial"), keep_prob=0.5, lambda=0, init=numeric(ncol(x)), gamma0=1, decay=0.01, n_passes=1000, minibatch_size=nrow(x))
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
  
  return(droplassoC(x, y, family, keep_prob, lambda, init, gamma0, decay, n_passes, minibatch_size))
}
