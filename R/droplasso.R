#' Fit a droplasso model
#' 
#' Fit a dropout lasso (droplasso) model. The regularization path is computed 
#' for the lasso component of the penalty at a grid of values for the 
#' regularization parameter lambda.
#' 
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is an 
#'   observation vector.
#' @param y Response variable, a vector of length \code{nobs} of quantitative 
#'   values for \code{family="gaussian"}, or of factors with two levels for 
#'   \code{family="binomial"}.
#' @param family Response type. \code{family="gaussian"} (default) for least 
#'   squares regression, \code{family="binomial"} for logistic regression
#' @param keep_prob The probability that each element is kept (default: 
#'   \code{0.5})
#' @param nlambda The number of \code{lambda} values (default: \code{10})
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
#' @param decay Learning rate decay (default: \code{1})
#' @param n_passes Number of passes over each example of the data on average 
#'   (default: \code{1000})
#' @param minibatch_size Batch size (default: \code{nobs})
#' @return An object of class \code{"droplasso"}, i.e. a list with the 
#'   following: \item{beta}{The \code{nvars x nlambda} matrix of weights, one 
#'   column per lambda, one row per variable} \item{lambda}{The sequence of 
#'   lambda for which the weigth is given} \item{nzero}{The number of non-zero 
#'   coefficient in each model} \item{call}{The function call}
#' @details Droplasso estimates a linear model by minimizing an objective 
#'   function \deqn{\min_{w} R(w) + \lambda*||w||_1} where \eqn{R(w)} is the 
#'   expected loss when the linear model is applied to a random training example
#'   subject to dropout noise, i.e., each coordinate is kept intact with 
#'   probability \code{keep_prob} and set to zero with probability \code{1 - 
#'   keep_prob}.
#' @details Given a prediction \eqn{u} and a true label \eqn{y}, the loss is
#'   \eqn{(u-y)^2 / 2} when \code{family="gaussian"}, and \eqn{-y*u + ln( 1+e^u
#'   )} when \code{family="binomial"} (i.e., the negative log-likelihood of the
#'   logistic regression model).
#' @details The optimization problem is solved with a stochastic proximal
#'   gradient descent algorithm, using mini-batches of size
#'   \code{minibatch_size}, and a learning rate decaying as 
#'   \code{gamma0/(1+decay*t)}, where \code{t} is the number of mini-batches 
#'   processed.
#' @details The problem is solved for all regularization parameters provided in
#'   the \code{lambda} argument. If no \code{lambda} argument is provided, then
#'   the function automatically chooses a decreasing sequence of \eqn{\lambda}'s
#'   to start from the null model and add features in the model along the
#'   regularization path. We use warm restart to start optimization for a given
#'   \eqn{\lambda} from the solution of the previous \eqn{lambda}, therefore it
#'   is strongly recommended to provide a sequence of \eqn{\lambda$}'s in
#'   decreasing order if a sequence is provided in the \code{lambda} argument.
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
#' @useDynLib droplasso
#' @importFrom Rcpp sourceCpp
droplasso <- function(x, y, family=c("gaussian","binomial"), keep_prob=0.5, nlambda=10, lambda.min.ratio = ifelse(nobs<nvars,0.01,0.0001), lambda=NULL, init=matrix(0,nrow=ncol(x)), gamma0=1, decay=1, n_passes=1000, minibatch_size=nrow(x))
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
  
  # Number of samples and features
  np=dim(x)
  if(is.null(np)|(np[2]<=1))stop("x should be a matrix with 2 or more columns")
  nobs=as.integer(np[1])
  nvars=as.integer(np[2])
  
  # Lambda values
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
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is an 
#'   observation vector.
#' @param y Response variable, a vector of length \code{nobs} of quantitative 
#'   values for \code{family="gaussian"}, or of factors with two levels for 
#'   \code{family="binomial"}.
#' @param family Response type. \code{family="gaussian"} (default) for least 
#'   squares regression, \code{family="binomial"} for logistic regression
lambda_max <- function(x, y, family=c("gaussian","binomial")) {
  family = match.arg(family)
  offset <- switch(family,
                   gaussian=0,
                   binomial=0.5)
  n <- dim(x)[1]
  s <- max(abs(t(x) %*% (y-offset)))/n
  return(s)
}
