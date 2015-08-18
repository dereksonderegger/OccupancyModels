#' Create Data using Royle-Nichols Occupancy Model
#' 
#' Blah blah blah
#' 
#' @param Design Data specifying the design of the data. A data frame with columns Species, Site, Day, Observer 
#' @param Abundance.formula A fixed effects formula for Abundance at each site
#' @param Abundance.params  Parameters that define Abundance
#' @param Detection.formula A fixed effects formula for Detection 
#' @param Detection.params  Parameters that define Detection
#' @return Somthing 
#' @export              
makeData.RN <- function(Design, Covariate.data=NULL,
                        Abundance.formula=~1,  Abundance.params=0,
                        Detection.formula=~1,  Detection.params=0){
  if(is.null(Covariate.data)){
    Covariate.data = data.frame(1:nrow(Design))
  }
  colnames(Design) <- c('Species','Site','Day','Observer')

  X <- model.matrix(Abundance.formula, Covariate.data)
  W <- model.matrix(Detection.formula, Covariate.data)
  out <- 
    Design %>% 
    mutate(
      lambda = as.vector(exp(       X %*% Abundance.params )),
      r      = as.vector(inv.logit( W %*% Detection.params )) ) %>%
    group_by(Species, Site) %>%
    mutate( N = rep( rpois(1, lambda[1]), n() )) %>%
    group_by(Species, Site, Day, Observer) %>%
    mutate( Detection = rbinom(1, 1, prob=1 - (1-r)^N) ) %>%
    arrange(Species, Site, Day, Observer)
  return(out)
} 
 


#' Royle-Nichols Occupancy Model
#' 
#' Blah blah some description of the model
#' 
#' @param Y  A data frame with 5 columns - Species, Site, Observer, Effort, NumDetections 
#' @param Covariate.data Data frame with the same number of rows as Y
#' @param Abundance.formula Model formula specifying the relationship between Abundance and the covariates
#' @param Detection.formula Model formula specifying the relationship between Detection probability and the covariates
#' @param N.limit Maximum conceivable number of individuals at a site. Making this too large slows down the model fitting, though.
#' @param n.chains How many MCMC chains to fit.
#' @param n.adapt How long should the burn-in phase be.
#' @param n.iter How many samples to take after the burn-in.
#' @param num.cores How many computer cores to use.
#' @param inits Initial values for the model chains.
#' @return An object of type RN_OccupancyModel
#' @export
Occ.RN <- function( Y, Covariate.data=NULL, 
                    Abundance.formula = ~1,
                    Detection.formula = ~1,
                    N.limit = 10,
                    n.chains=4, n.adapt=1000, n.iter=1000,
                    num.cores=1,
                    fit=NA, inits=NULL ){
  
  colnames(Y) <- c('..Species','..Site','..Camera','..Effort','..NumDetections')
  Y <- Y %>% group_by() %>%
    mutate(..Species = as.factor(..Species), 
           ..Site    = as.factor(..Site),
           ..Camera  = as.factor(..Camera) ) %>%
    droplevels() %>%
    mutate(..SpeciesSite = interaction(..Species, ..Site)) %>%
    mutate_each(funs(as.integer)) %>%
    select(..Species, ..Site, ..Camera, ..Effort, ..NumDetections )
  
  data <- cbind(Y, Covariate.data)
  data <- droplevels(data)  # drop unused factor levels
  
  Data <- list(Y=Y, n=nrow(Y))
  pars <- NULL
  if( is.null( findbars(Abundance.formula) ) ){
    Data$AX   <- model.matrix( Abundance.formula, data )
    Data$p_AX <- ncol(Data$AX)
    Data$AZ   <- model.matrix( ~1, data)
    Data$p_AZ <- 1
    Data$q_AZ <- 1
    Data$i_AZ_size <- 2
    Data$i_AZ <- c(1,1)  # can be anything but stan wants it to be a vector
    pars <- c(pars, 'beta_A')
  }else{
    temp    <- lFormula( as.formula( paste('..NumDetections', deparse(Abundance.formula))), data )
    Data$AX   <- temp$X
    Data$p_AX <- ncol(Data$AX)
    Data$AZ   <- as.matrix( t(temp$reTrms$Zt) )
    Data$p_AZ <- ncol(Data$AZ)
    Data$q_AZ <- max(temp$reTrms$Lind)
    Data$i_AZ <- temp$reTrms$Lind
    Data$i_AZ_size <- Data$p_AZ
    pars <- c(pars, 'beta_A', 'gamma_A', 'sigma_A')
  }

  if( is.null(findbars(Detection.formula)) ){
    Data$DX   <- model.matrix( Detection.formula, data )
    Data$p_DX <- ncol(Data$DX)
    Data$DZ   <- model.matrix( ~1, data)
    Data$p_DZ <- 1
    Data$q_DZ <- 1
    Data$i_DZ_size <- 2
    Data$i_DZ <- c(1,1)  # can be anything but stan wants it to be a vector
    pars <- c(pars, 'beta_D')
  }else{
    temp    <- lFormula( as.formula( paste('..NumDetections', deparse(Detection.formula))), data )
    Data$DX   <- temp$X
    Data$p_DX <- ncol(Data$DX)
    Data$DZ   <- as.matrix( t(temp$reTrms$Zt) )
    Data$p_DZ <- ncol(Data$DZ)
    Data$q_DZ <- max(temp$reTrms$Lind)
    Data$i_DZ <- temp$reTrms$Lind
    Data$i_DZ_size <- Data$p_DZ
    pars <- c(pars, 'beta_D', 'gamma_D', 'sigma_D')
  }
  #pars <- c(pars, 'r')
  Data$N_limit <- N.limit
  
  Occ_compiled_RN_Model_STAN <<- stan_model( file=system.file('RN_Model.stan', package='OccupancyModels'),
                                             model_name='Occ_RN_Model_STAN' )
  chains <- sampling(Occ_compiled_RN_Model_STAN, data=Data,
                  pars=pars, chains=n.chains, iter=n.iter, cores=num.cores )
#   out <- stan( model_code=RN_Model(), pars = pars, data = Data, 
#                chains = n.chains, iter=n.iter)

  Data$Abundance.formula <- Abundance.formula
  Data$Detection.formula <- Detection.formula
  Data$Covariate.data <- Covariate.data
  
    
  # If N.limit is less than, say 10*max lambda, we should
  # throw a warning and tell the user to re-run the model
  # with a higher N.limit
#   lambda.random <- extract(out, pars='lambda_random')$lambda_random
#   max.lambda <- max(lambda.random)
#   min.N.limit <- qpois(.99, max.lambda)
#   if( N.limit < min.N.limit){
#     warning(str_c('The maximum expected abuandance was ', round(max.lambda, digits=1), 
#                   ' and N.limit should be at least ', round(min.N.limit)) )
#   }
  
  out2 <- as.RN_OccupancyModel(Data, 
                               stan.model = Occ_compiled_RN_Model_STAN,
                               chains = chains,
                               pars)
  

  return(out2)  
}




#' Royle-Nichols Model convert beta to psi
#' 
#' @param beta A vector of length p
#' @param X A design matrix of dimensions n-by-p
#' @export
RN.beta2Psi <- function(beta, X){
  return( 1 - exp(-exp(X %*% beta)) )
}

#' Royle-Nichols Model convert psi to beta
#' 
#' @param Psi A vector of length n
#' @param X A design matrix of dimensions n-by-p
#' @export
RN.Psi2beta <- function(Psi, X){
  return( solve(t(X)%*%X) %*% t(X)%*%log(-log(1-Psi)) )
}

#' Royle-Nichols Model convert alpha to r
#' 
#' @param alpha A vector of length p
#' @param W A design matrix of dimensions n-by-p
#' @export
RN.alpha2r  <- function(alpha, W){
  return( inv.logit(W %*% alpha))
}

#' Royle-Nichols Model convert r to alpha
#' 
#' @param r A vector of length n
#' @param W A design matrix of dimensions n-by-p
#' @export
RN.r2alpha  <- function(r, W){
  return(solve(t(W)%*%W) %*% t(W) %*% logit(r))
}

