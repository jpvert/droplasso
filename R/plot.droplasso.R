#' Plot coefficients from a "droplasso" object
#' 
#' @param m Fitted model of class "survenet"
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
plot.droplasso = function(m, sign.lambda = 1, ...) {
  loglambda = log(m$lambda)
  xlab = "Log lambda"
  if (sign.lambda < 0)
    xlab = paste("-", xlab, sep = "")
  index = sign.lambda * loglambda
  matplot(index,t(m$beta),lty=1,xlab=xlab,ylab="Coefficients",type="l",...)
}