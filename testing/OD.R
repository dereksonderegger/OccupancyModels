# A simple case making sure we can run the model
n <- 100;   Psi <- .8;   r <- .8;   Effort <- 15
# make some data
Design <- expand.grid(Species=1, Site=factor(1:n), Day=factor(1:Effort), Observer='C1')
Covariate.data <- data.frame( group=rep(c('A'), each=n))
X <- model.matrix(~ 1, Covariate.data)
W <- model.matrix(~ 1, Covariate.data)
beta_OX <- OD.Psi2beta(rep(Psi,nrow(X)), X)
beta_DX <- OD.r2beta( rep(r,  nrow(W)), W)

data_OD_long <- makeData.OD(
  Design, 
  Occupancy.formula = ~ 1,  Occupancy.params = beta_OX,
  Detection.formula = ~ 1,  Detection.params = beta_DX  )
data_OD_short <- data_OD_long %>% group_by(Species, Site, Observer) %>% 
  summarize( Effort=n(), Detections=sum(Detection) ) %>%
  select(Species, Site, Observer, Effort, Detections)
# Sanity check if the data is good
data_OD_short %>% group_by(Species, Site) %>% 
  summarize(Occupied = as.integer(max(Detections)>0)) %>%
  summarize(Percent.Occupied = sum(Occupied)/n())

object <- Occ.OD(   Y = data_OD_short, 
                    Covariate.data = Covariate.data,
                    Occupancy.formula = ~ 1,
                    Detection.formula = ~ 1,
                    n.chains=4, n.iter=1000)

summary(object)
beta_OX; beta_DX
traceplot(object$chains, pars='beta_D')

ggs_traceplot(object$chains, family='beta')
object <- update.Occ.chain(object, n.iter=1000)
ggs_traceplot(object$chains, family='beta')

summary(object)



############################
## A multi-species test
############################
n <- 300;   Effort <- 50
# make some data
Design <- expand.grid(
  Species=c('Coyote','Deer','Cougar'),
  Site=factor(1:n), 
  Day=factor(1:Effort), 
  # Observer=c('OnTrail','OffTrail') )
  Observer=c('OnTrail') ) %>%
  arrange(Species, Site, Day, Observer)
  
Site.Data <- data.frame(
  Site    = factor(1:n),
  Habitat = sample( c('Forest','Plains','Desert'), n, replace=TRUE) )

covariates <- merge(Design, Site.Data) %>% 
  arrange(Species, Site, Day, Observer)

# X <- model.matrix(~ Species+Habitat, covariates)
X <- model.matrix(~ Species, covariates)
# W <- model.matrix(~ Observer, covariates)
W <- model.matrix(~ 1, covariates)

#               Coyote,Desert  +Deer   +Cougar  +Forest  +Plains
#Psi <- X %*% c(      .4,        .2,     -.35,      .2,      .1)
Psi <- X %*% c(      .4,        .4,     -.2)
beta <- OD.Psi2beta(Psi, X)

#                OnTrail  +OffTrail
#r   <- W %*% c(  .8,        -.6 )
r   <- W %*% c(  .8 )
alpha <- OD.r2beta(r , W)

data_OD_long <- makeData.OD(
  Design, 
  Covariate.data=covariates,
  Occupancy.formula = ~ Species,  
  Detection.formula = ~ 1,  
  Occupancy.params = beta,
  Detection.params = alpha) 
data_OD_short <- data_OD_long %>% group_by(Species, Site, Observer) %>% 
  summarize( Effort=n(), Detections=sum(Detection) ) %>%
  arrange(Species, Site, Observer) %>%
  select(Species, Site, Observer, Effort, Detections)
data <- merge(data_OD_short, Site.Data) %>%
  arrange(Species, Site, Observer)

data %>% group_by(Species, Site) %>% 
  summarize(Occupied = as.integer(max(Detections)>0)) %>%
  summarize(Percent.Occupied = sum(Occupied)/n())


object <- Occ.OD( Y = data %>% select(Species,Site,Observer,Effort,Detections),
                  Covariate.data = data,
                  Occupancy.formula = ~ Species,
                  Detection.formula = ~ 1,
                  n.chains=4, n.iter=1000)

summary(object)
rstan::traceplot(object$chains, pars='beta_O')
beta; alpha

