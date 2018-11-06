#' Make predictions from a droplasso object.
#' 
#' Similar to other \code{predict} methods, this function predicts fitted values
#' or coefficients from a fitted "\code{droplasso}" object
#' @param object Fitted "\code{droplasso}" model object
#' @param newx Matrix of new values for \code{x} at which predictions are to be 
#'   made. This argument is not used for \code{type="coefficients"}.
#' @param s Value(s) of the penalty parameter lambda at which predictions are 
#'   required. Can be either omitted, or a vector of numeric values. If omitted
#'   (or equal to \code{NULL}), default is the entire sequence used to create
#'   the model.
#' @param type Type of prediction required. Type \code{"link"} (default) gives 
#'   the linear predictor. Type \code{"coefficients"} computes the coefficients 
#'   at the requested values of \code{s}. If \code{newx} is not provided, then 
#'   \code{type="coefficients"} automatically.
#' @return If \code{type="link"}, returns a matrix of prediction scores, one row
#'   per sample of \code{newx}, one column per value of lambda. If 
#'   \code{type="coefficients"}, returns a matrix of weights of the model, one 
#'   row per feature, one column per value of lambda.
#' @seealso \code{droplasso}, \code{coef} methods.
#' @examples
#' #create data:
#' nobs = 100
#' nvars = 5
#' x = matrix(rnorm(nobs*nvars),nrow=nobs)
#' b = c(1,1,0,0,0)
#' p = 1/(1+exp(-x%*%b))
#' y = p>0.5
#' # Fit a dropout lasso model
#' m <- droplasso(x, y, family="binomial", nlambda=50, decay=5)
#' # Plot is
#' plot(m)
#' # Plot the weights for lambda=0.01
#' plot(predict(m,s=0.01), xlab="Features", ylab="Weight", main="Lambda=0.01")
#' # Plot prediction score
#' plot(predict(m, x, s=0.01), xlab="Sample", ylab="Score", main="Lambda=0.01")
#' @export
predict.droplasso = function(object, newx, s=NULL, type=c("link", "coefficients")) {
  
  type = match.arg(type)
  if (missing(newx)) {
    type="coefficients"
  }
  
  if (!is.null(s)) {
    if (!is.numeric(s)){stop("Invalid form for s")}
  }
  
  nbeta=object$beta
  if (!is.null(s)) {
    # If s is not the default, make the matrix of predictors for the corresponding values of lambda
    vnames = dimnames(nbeta)[[1]]
    dimnames(nbeta) = list(NULL, NULL)
    lambda = object$lambda
    lamlist = lambda.interp(lambda, s)
    nbeta = nbeta[, lamlist$left, drop = FALSE] %*% diag(lamlist$frac,length(s)) + nbeta[, lamlist$right, drop = FALSE] %*% diag(1 - lamlist$frac,length(s))
    dimnames(nbeta) = list(vnames, paste(seq(along = s)))
  }
  if (type == "coefficients")
    return(nbeta)

  nfit = as.matrix(newx %*% nbeta)
  return(nfit)
}



# Interpolate s values on a sequence of lambda
# Taken from glmnet package
lambda.interp = function(lambda, s)
{
  if (length(lambda) == 1) {
    nums = length(s)
    left = rep(1, nums)
    right = left
    sfrac = rep(1, nums)
  }
  else {
    s[s > max(lambda)] = max(lambda)
    s[s < min(lambda)] = min(lambda)
    k = length(lambda)
    sfrac <- (lambda[1] - s)/(lambda[1] - lambda[k])
    lambda <- (lambda[1] - lambda)/(lambda[1] - lambda[k])
    coord <- approx(lambda, seq(lambda), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac = (sfrac - lambda[right])/(lambda[left] - lambda[right])
    sfrac[left == right] = 1
  }
  list(left = left, right = right, frac = sfrac)
}
