#' Create simulated Occupancy data using the Occupancy/Detection model 
#' 
#' Creates simulated data from the standard
#' Occupancy/Availability model. 
#' @param Design Data specifying the design of the data. A data frame with columns Species, Site, Day, Observer 
#' @param Covariate.data 
#' @param Occupancy model A fixed effects formula for Occupancy at each site
#' @param Occupancy params Parameters for the Occupancy formula
#' @param Detection model A fiexed effects formula for Detection at each site
#' @param Detection params Parameters for the Detection formula
#' @export
makeData.OD <- function(Design, Covariate.data = NULL,
                      Occupancy.formula = ~1, Occupancy.params=c(0),
                      Detection.formula = ~1, Detection.params=c(0) ){

  if(is.null(Covariate.data)){
    Covariate.data = data.frame(1:nrow(Design))
  }
  colnames(Design) <- c('Species','Site','Day','Observer')
  
  OX <- model.matrix(Occupancy.formula, Covariate.data)
  DX <- model.matrix(Detection.formula, Covariate.data)
  
  out <- Design %>% mutate(
            Psi = as.vector(inv.logit(OX %*% Occupancy.params)),
            r = as.vector(inv.logit(DX %*% Detection.params))) %>%
         group_by(Species,Site) %>%
         mutate( Occupied = rep( rbinom(1,1,prob=Psi), n() )) %>%
         group_by(Site, Day, Observer) %>%
         mutate( Detection = Occupied * rbinom(n(), 1, prob=r) )

  return(out) 
}


#' Analyzes Occupancy data using the Occupancy/Detection model 
#' 
#' Uses Mackenzie's Occupancy/Detection model.
#' @param Y  A data frame with 5 columns - Species, Site, Observer, Effort, NumDetections 
#' @param Covariate.data Data frame with the same number of rows as Y
#' @param Occupancy.formula Model formula specifying the relationship between Occupancy and the covariates
#' @param Detection.formula Model formula specifying the relationship between Detection probability and the covariates
#' @param N.limit Maximum conceivable number of individuals at a site. Making this too large slows down the model fitting, though.
#' @param n.chains How many MCMC chains to fit.
#' @param n.adapt How long should the burn-in phase be.
#' @param n.iter How many samples to take after the burn-in.
#' @param num.cores How many computer cores to use.
#' @param inits Initial values for the model chains.
#' @param ... Further arguments passed to Stan's sampling() command
#' @export
Occ.OD <- function(Y, Covariate.data=NULL, 
                   Occupancy.formula = ~1,
                   Detection.formula = ~1,
                   n.chains=4, n.adapt=1000, n.iter=1000,
                   num.cores=1,
                   fit=NULL, inits=NULL, ...){

  colnames(Y) <- c('..Species','..Site','..Camera','..Effort','..NumDetections')
  Y <- Y %>% group_by() %>%
    mutate(..Species = as.factor(..Species), 
           ..Site    = as.factor(..Site),
           ..Camera  = as.factor(..Camera),
           ..KnownOccupied = as.integer(..NumDetections >= 1)) %>%
    droplevels() %>%
    mutate(..SpeciesSite = interaction(..Species, ..Site)) %>%
    mutate_each(funs(as.integer)) %>%
    select(..Species, ..Site, ..Camera, ..KnownOccupied, ..Effort, ..NumDetections )
  
  
  data <- cbind(Y, Covariate.data)
  data <- droplevels(data)  # drop unused factor levels
  
  Data <- list(Y=Y, n=nrow(Y))
  pars <- NULL
  
  if( is.null( findbars(Occupancy.formula) ) ){
    Data$OX   <- model.matrix( Occupancy.formula, Covariate.data )
    Data$p_OX <- ncol(Data$OX)
    Data$OZ   <- model.matrix( ~1, Covariate.data)
    Data$p_OZ <- 1
    Data$q_OZ <- 1
    Data$i_OZ_size <- 2
    Data$i_OZ <- c(1,1)  # can be anything but stan wants it to be a vector
    pars <- c(pars, 'beta_O')
  }else{
    temp    <- lFormula( update.formula(Occupancy.formula, ..Occ.Y~.), Covariate.data )
    Data$OX   <- temp$X
    Data$p_OX <- ncol(Data$OX)
    Data$OZ   <- as.matrix( t(temp$reTrms$Zt) )
    Data$p_OZ <- ncol(Data$OZ)
    Data$q_OZ <- max(temp$reTrms$Lind)
    Data$i_OZ <- temp$reTrms$Lind
    Data$i_OZ_size <- Data$p_OZ
    pars <- c(pars, 'beta_O', 'gamma_O', 'sigma_O')
  }
  if( is.null( findbars(Detection.formula) )){
    Data$DX   <- model.matrix( Detection.formula, Covariate.data )
    Data$p_DX <- ncol(Data$DX)
    Data$DZ   <- model.matrix( ~1, Covariate.data)
    Data$p_DZ <- 1
    Data$q_DZ <- 1
    Data$i_DZ_size <- 2
    Data$i_DZ <- c(1,1)  # can be anything but stan wants it to be a vector
    pars <- c(pars, 'beta_D')
  }else{
    temp    <- lFormula( update.formula(Detection.formula, ..Occ.Y~.), Covariate.data )
    Data$DX   <- temp$X
    Data$p_DX <- ncol(Data$DX)
    Data$DZ   <- as.matrix( t(temp$reTrms$Zt) )
    Data$p_DZ <- ncol(Data$DZ)
    Data$q_DZ <- max(temp$reTrms$Lind)
    Data$i_DZ <- temp$reTrms$Lind
    Data$i_DZ_size <- Data$p_DZ
    pars <- c(pars, 'beta_D', 'gamma_D', 'sigma_D')
  }

  Occ_compiled_OD_model_STAN <<- 
    stan_model( file=system.file('OD_Model.stan', package='OccupancyModels'),
                model_name = 'Occ_OD_STAN')
  
  chains <- sampling(Occ_compiled_OD_model_STAN, data=Data,
            pars=pars, chains=n.chains, iter=n.iter, cores=num.cores, ... )
  
  Data$Occupancy.formula <- Occupancy.formula
  Data$Detection.formula <- Detection.formula
  Data$Covariate.data <- Covariate.data
  
  out2 <- as.OD_OccupancyModel(Data, 
                               stan.model = Occ_compiled_OD_model_STAN,
                               chains = chains,
                               pars)

  return(out2)
}


#' Mckenzie Model convert beta to psi
#' 
#' @param beta A vector of length p
#' @param X A design matrix of dimensions n-by-p
#' @export
OD.beta2Psi <- function(beta, X){
  return( inv.logit(X %*% beta))
}

#' Mckenzie Model convert psi to beta
#' 
#' @param Psi A vector of length n
#' @param X A design matrix of dimensions n-by-p
#' @export
OD.Psi2beta <- function(Psi, X){
  return(solve(t(X)%*%X) %*% t(X) %*% logit(Psi))
}

#' McKenzie Model convert alpha to r
#' 
#' @param alpha A vector of length p
#' @param W A design matrix of dimensions n-by-p
#' @export
OD.beta2r  <- function(beta, W){
  return( inv.logit(W %*% beta))
}

#' McKenzie Model convert r to alpha
#' 
#' @param r A vector of length n
#' @param W A design matrix of dimensions n-by-p
#' @export
OD.r2beta  <- function(r, W){
  return(solve(t(W)%*%W) %*% t(W) %*% logit(r))
}


