data {
  int <lower=1> n;   // Number of sites/observations
  int Y[n,6];        // Species, SiteID, CameraID, SiteKnownOccupied, Effort, NumDetections
  int p_OX;          // Number of X columns (Occupancy)
  int p_OZ;          // Number of Z columns random effects (Occupancy)
  int q_OZ;          // Number of sigma values (Occupancy)
  int p_DX;          // Number of X columns (Detection)
  int p_DZ;          // Number of z columns random effects (Detection)
  int q_DZ;          // Number of sigma values (Detection)
  matrix [n,p_OX] OX;   // Occupancy fixed design matrix
  matrix [n,p_OZ] OZ;   // Occupancy Random design matrix
  matrix [n,p_DX] DX;   // Detection fixed design Matrix
  matrix [n,p_DZ] DZ;   // Detection Random design matrix
  int i_OZ_size;
  int i_DZ_size;
  int i_OZ[i_OZ_size];  // Abundance sigma -> gamma design matrix
  int i_DZ[i_DZ_size];  // Detection sigma -> gamma design matrix
}
parameters {
  vector <lower=-5, upper=5> [p_OX] beta_O;
  vector <lower=-5, upper=5> [p_OZ] gamma_O;
  vector <lower=-5, upper=5> [p_DX] beta_D;
  vector <lower=-5, upper=5> [p_DZ] gamma_D;
  vector <lower=0, upper=10> [q_OZ] sigma_O;
  vector <lower=0, upper=10> [q_DZ] sigma_D;
}
transformed parameters{
  vector [n] Psi;
  vector [n] r;
  vector [n] logit_Psi;
  vector [n] logit_r;
  if( p_OZ > 1 ){
    logit_Psi  <-  OX * beta_O + OZ * gamma_O;
  }else{
    logit_Psi  <- OX * beta_O;
  }
  if( p_DZ > 1){
    logit_r <- DX * beta_D + DZ * gamma_D;
  }else{
    logit_r <- DX * beta_D;
  }
  for(i in 1:n){
    Psi[i] <- inv_logit(logit_Psi[i]);
    r[i]   <- inv_logit(logit_r[i]);
  }
}
model{
  real loglikely_cameras;
  int  next_site;
  int  next_day;

   beta_O ~ normal(0,10);
   beta_D ~ normal(0,10);
  sigma_O ~ cauchy(0,5);
  sigma_D ~ cauchy(0,5);
  gamma_O ~ normal(0, sigma_O);
  gamma_D ~ normal(0, sigma_D);
      
  loglikely_cameras <- 0;
  for( i in 1:n ){
    if( i == n ){
      next_site <- Y[i,2]+1;
    }else{
      next_site <- Y[i+1,2];
    }
      
    // update camera likelyhood      
    loglikely_cameras <- loglikely_cameras + binomial_log(Y[i,6],Y[i,5],r[i]);
    #likely_cameras * exp(binomial_coefficient_log(Y[i,5],Y[i,6])) *
    #                r[i]^(Y[i,4]) * (1-r[i])^(Y[i,5]-Y[i,6]); 
    if( Y[i,1] != next_site ){
      increment_log_prob( 
        log(  (1-Psi[i])*(1-Y[i,4]) + Psi[i]*exp(loglikely_cameras) ));
      loglikely_cameras <- 0;
    }
  }
}
