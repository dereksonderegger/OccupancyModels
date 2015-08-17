
#' Create a RN_OccupancyModel object
as.RN_OccupancyModel <- function(data, stan.model, chains, pars){
  out <- list(user.data = data,
              stan.model = stan.model,
              chains = chains,
              pars = pars)
  class(out) <- 'RN_OccupancyModel'
  return(out)
}

#' @export
traceplot.RN_OccupancyModel <- function(object, pars=NULL){
  traceplot(object$chains, pars)
}

#' Return the model parameter estimates
coef.RN_OccupancyModel <- function(obj, family=NULL){
  output <- get_family(obj, family) %>%
    group_by(Parameter) %>%
    summarise( Est = mean(value))
  return(as.data.frame(output))
}


#' @export
predict.RN_OccupancyModel <- function(object, newdata=NULL, level=0.95){
  # Give the Psi values along with upper and lower CI
  if(is.null(newdata)){
    AX <- object$user.data$AX
  }else{
    AX <- model.matrix(object$user.data$Occupancy.formula, newdata) 
  }
  lambda <- exp(AX %*% t(rstan::extract(object$chains, pars='beta_A')$beta_A) )
  lambda <- cbind( 
    lambda.hat = 1-ppois(0, apply(lambda, 1, mean)),
    lwr        = 1-ppois(0, apply(lambda, 1, quantile, probs=   (1-level)/2)),
    upr        = 1-ppois(0, apply(lambda, 1, quantile, probs= 1-(1-level)/2)) )
  
  return(Psi)    
}


#' @export
confint.RN_OccupancyModel <- function(object, fun=identity, level=.95){
  temp <- summary.RN_OccupancyModel(object, fun=fun, level=level)
  output <- NULL
  output$Abundance <- temp$Abundance %>%
    filter(Parameter, lwr, upr)
  output$Detection <- temp$Detection %>%
    filter(Parameter, lwr, upr)
  return(output)
}

#' @export
coef.RN_OccupancyModel <- function(object, fun=identity){
  temp <- summary.RN_OccupancyModel(object, fun=fun, level=level)
  output <- NULL
  output$Abundance <- temp$Abundance %>%
    filter(Parameter, median)
  output$Detection <- temp$Detection %>%
    filter(Parameter, median)
  return(output)
}


#' @export
summary.RN_OccupancyModel <- function(object, fun=identity, level=0.95){
  
  alpha <- 1-level
  tails <- c(alpha/2, 1-alpha/2)
  composite<-function(f,g) function(...) f(g(...))
  
  output <- NULL
  
  # Figure out the Abundance Parameters
  betas <- as.data.frame(rstan::extract(object$chains, pars='beta_A')$beta_A) 
  colnames(betas) <- colnames(object$user.data$AX)
  betas <- gather(as.data.frame(betas), key='Parameter')
  output$Abundance <- betas %>% group_by(Parameter) %>%
  #betas <- get_family(object, family='beta_A')  
  #output$Occupancy <- betas %>% group_by(Parameter) %>%
    mutate( value = fun(value)) %>%
    summarise( mean = mean(value),
               median = median(value),
               lwr = quantile(value, probs=alpha/2),
               upr = quantile(value, probs=1-alpha/2))
  
  # Figure out the Detection Parameters
  betas <- as.data.frame(rstan::extract(object$chains, pars='beta_D')$beta_D) 
  colnames(betas) <- colnames(object$user.data$DX)
  betas <- gather(as.data.frame(betas), key='Parameter')
  ##betas <- get_family(object, family='beta_D')  
  output$Detection <- betas %>% group_by(Parameter) %>%
    mutate( value = fun(value)) %>%
    summarise( mean = mean(value),
               median = median(value),
               lwr = quantile(value, probs=alpha/2),
               upr = quantile(value, probs=1-alpha/2)) 
  
  return(output)
}

