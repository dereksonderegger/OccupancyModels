
#' Create a OD_OccupancyModel object
as.OD_OccupancyModel <- function(data, stan.model, chains, pars){
  out <- list(user.data = data,
              stan.model = stan.model,
              chains = chains,
              pars = pars)
  class(out) <- 'OD_OccupancyModel'
  return(out)
}

#' @export
traceplot.OD_OccupancyModel <- function(object, pars=NULL){
  traceplot(object$chains, pars)
}

#' Return the model parameter estimates
coef.OD_OccupancyModel <- function(obj, family=NULL){
  output <- get_family(obj, family) %>%
    group_by(Parameter) %>%
    summarise( Est = mean(value))
  return(as.data.frame(output))
}


#' @export
predict.OD_OccupancyModel <- function(object, newdata=NULL, level=0.95){
  # Give the Psi values along with upper and lower CI
  if(is.null(newdata)){
    OX <- object$user.data$OX
  }else{
    OX <- model.matrix(object$user.data$Occupancy.formula, newdata) 
  }
  Psi <- OX %*% t(rstan::extract(object$chains, pars='beta_O')$beta_O) 
  Psi <- inv.logit(cbind( 
    Psi.hat = apply(Psi, 1, mean),
    lwr     = apply(Psi, 1, quantile, probs=   (1-level)/2),
    upr     = apply(Psi, 1, quantile, probs= 1-(1-level)/2) ))

  return(Psi)    
}

#' @export
confint.OD_OccupancyModel <- function(object, fun=identity, level=.95){
  temp <- summary.OD_OccupancyModel(object, fun=fun, level=level)
  output <- NULL
  output$Occupancy <- temp$Occupancy %>%
    filter(Parameter, lwr, upr)
  output$Detection <- temp$Detection %>%
    filter(Parameter, lwr, upr)
  return(output)
}

#' @export
coef.OD_OccupancyModel <- function(object, fun=identity){
  temp <- summary.OD_OccupancyModel(object, fun=fun, level=level)
  output <- NULL
  output$Occupancy <- temp$Occupancy %>%
    filter(Parameter, median)
  output$Detection <- temp$Detection %>%
    filter(Parameter, median)
  return(output)
}


#' @export
summary.OD_OccupancyModel <- function(object, fun=identity, level=0.95){

  alpha <- 1-level
  tails <- c(alpha/2, 1-alpha/2)
  composite<-function(f,g) function(...) f(g(...))
  
  output <- NULL
  
  # Figure out the Occupancy Parameters
  betas <- as.data.frame(rstan::extract(object$chains, pars='beta_O')$beta_O) 
  colnames(betas) <- colnames(object$user.data$OX)
  betas <- gather(as.data.frame(betas), key='Parameter')
  output$Occupancy <- betas %>% group_by(Parameter) %>%
    mutate( value = fun(value)) %>%
    summarise( mean = mean(value),
               median = median(value),
               lwr = quantile(value, probs=alpha/2),
               upr = quantile(value, probs=1-alpha/2)) 
  
  # Figure out the Detection Parameters
  betas <- as.data.frame(rstan::extract(object$chains, pars='beta_D')$beta_D) 
  colnames(betas) <- colnames(object$user.data$DX)
  betas <- gather(as.data.frame(betas), key='Parameter')
  output$Detection <- betas %>% group_by(Parameter) %>%
    mutate( value = fun(value)) %>%
    summarise( mean = mean(value),
               median = median(value),
               lwr = quantile(value, probs=alpha/2),
               upr = quantile(value, probs=1-alpha/2)) 
  
  return(output)
}

