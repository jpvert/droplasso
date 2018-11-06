#' Plot coefficients from a "droplasso" object
#' 
#' @param x Fitted model of class "survenet"
#' @param sign.lambda Should be \code{sign.lambda = 1} (default) to plot from
#'   left to right in increasing values of \code{lambda}, \code{sign.lambda =
#'   -1} otherwise
#' @param ... Other arguments passed to \code{\link[graphics]{matplot}}
#' @return None
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
#' @export
#' @importFrom graphics matplot
plot.droplasso = function(x, sign.lambda = 1, ...) {
  loglambda = log(x$lambda)
  xlab = "Log lambda"
  if (sign.lambda < 0)
    xlab = paste("-", xlab, sep = "")
  index = sign.lambda * loglambda
  matplot(index,t(x$beta),lty=1,xlab=xlab,ylab="Coefficients",type="l",...)
}