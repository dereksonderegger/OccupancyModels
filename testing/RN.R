n <- 40;   Psi <- .8;   r <- .8;   Effort <- 14
# make some data
Design <- expand.grid(Species=1, Site=factor(1:n), Day=factor(1:Effort), Observer='C1')
Covariate.data <- data.frame( group=rep(c('A'), each=n))
X <- model.matrix(~ 1, Covariate.data)
W <- model.matrix(~ 1, Covariate.data)
beta_A <- RN.Psi2beta(rep(Psi,nrow(X)), X)
beta_D <- RN.r2alpha(rep(r,nrow(X)), W)
data_RN_long <- makeData.RN(
  Design, 
  Abundance.formula = ~ 1,  Abundance.params = beta_A,
  Detection.formula = ~ 1,  Detection.params = beta_D )
data_RN_short <- data_RN_long %>% group_by(Species, Site, Observer, N) %>% 
  summarize( Effort=n(), Detections=sum(Detection) ) %>%
  select(Species, Site, Observer, Effort, Detections)
# Sanity check if the data is good
data_RN_short %>% group_by(Species, Site) %>% 
  summarize(Occupied = as.integer(max(Detections)>0)) %>%
  summarize(Percent.Occupied = sum(Occupied)/n())

data_RN_short$X <- rnorm(nrow(data_RN_short))

object <- Occ.RN(Y = select(data_RN_short, 
                            Species,Site,Observer,Effort,Detections),
                    Covariate.data = data_RN_short,
                    Abundance.formula = ~ X,
                    Detection.formula = ~ 1,
                    n.chains=4, n.iter=1000, num.cores=2)
summary(object)

rstan::traceplot(object$chains, pars='beta_D')
predict(object)
beta_A; beta_D






