#' Extend the MCMC chains
#' @export
update.Occ.chain <- function(object, n.iter){
  out <- coda.samples(object$model, object$pars, n.iter)
  object$chains <- ggs(out)
  return(object)
}

#' Return the mcmc chain as a tbl_df object
#' @export
get_family <- function(object, family=NULL){
  return( ggs_density(object$chains, family=family)$data )
}

#' Return the mcmc chain as a tbl_df object
#' @export
get.family <- function(object, family=NULL){
  return( ggs_density(object$chains, family=family)$data )
}


# Little function to calculate posterior variances from simulation
colVars <- function (a, na.rm=FALSE){
  diff <- a - matrix (colMeans(a, na.rm=na.rm), nrow(a), ncol(a), byrow=TRUE)
  vars <- colMeans (diff^2, na.rm=na.rm)*nrow(a)/(nrow(a)-1)
  return (vars)
}
waic <- function (stanfit, na.rm=FALSE){
  log_lik <- rstan::extract(stanfit, "log_lik")$log_lik
  lppd <- sum( log( colMeans( exp(log_lik), na.rm=na.rm)))
  p_waic_1 <- 2*sum (log(colMeans(exp(log_lik), na.rm=na.rm)) - colMeans(log_lik, na.rm=na.rm))
  p_waic_2 <- sum (colVars(log_lik, na.rm=na.rm))
  waic_2 <- -2*lppd + 2*p_waic_2
  return (list (waic=waic_2, p_waic=p_waic_2, lppd=lppd, p_waic_1=p_waic_1))
}




#' Convert Time from 0-24:00 to 0-2*pi
#' 
#' Because the package overlap parameterizes time of day 
#' from 0 to 2*pi, times need to be converted. 
#' @param x A POSIXct date/time object
#' @export
#' @examples 
#' x <- strptime('09-18-2010 17:30:05', format='%m-%d-%Y %H:%M:%S')
#' toRadianTime(x)
#' 
#' # You can also use package lubridate to create date/time objects
#' library(lubridate)
#' y <- mdy_hms('09-18-2010 06:00:05')
#' toRadianTime(y)
toRadianTime <- function(x){
  out <- as.numeric(format(x,'%H')) / ( 24 )     +
         as.numeric(format(x,'%M')) / ( 24*60 )  +
         as.numeric(format(x,'%S')) / ( 24*60*60)
  return(out * 2 * pi)
}



#' Logit function
#' 
#' Calculates the log odds ratio.
#' @param p a probability or proportion
#' @export
#' @examples
#' logit(.4)
logit <- function(p){
  log( p/(1-p) )
}

#' Inverse Logit function
#' 
#' @param x a real number 
#' @export
#' @examples
#' inv.logit(3)
inv.logit <- function(x){
  1/( 1+exp(-x))
}


#' Safe Averaging
#' 
#' A function that calculates the mean if given numerical values.
#' If passed a vector of factors, it checks if the levels are all the same
#' and returns that level, otherwise it can't average factor levels and
#' errors out.
#' @param x a vector of numerical, logical, or factors
#' @export
#' @examples
#' safe_mean( c(1,2,3) )
#' safe_mean( factor('A','A','A') ) 
#' # Not run, should throw an error.
#' # safe_mean( factor('A','B','B') )
safe_mean <- function(x){
  if( is.numeric(x) )
    return(mean(x))
  else if( length(unique(x)) == 1)
    return( x[1] )
  else
    stop(paste('Error safe_mean: Trying to average multiple levels.'))
}


