  data {
    int <lower=1> n;   // Number of sites/observations
    int Y[n,5];        // Species, SiteID, CameraID, Effort, NumDetections
    int p_AX;          // Number of X columns (Abundance)
    int p_AZ;          // Number of Z columns random effects (Abundance)
    int q_AZ;          // Number of sigma values (Abundance)
    int p_DX;          // Number of X columns (Detection)
    int p_DZ;          // Number of z columns random effects (Detection)
    int q_DZ;          // Number of sigma values (Detection)
    matrix [n,p_AX] AX;   // Abundance fixed design matrix
    matrix [n,p_AZ] AZ;   // Abundance Random design matrix
    matrix [n,p_DX] DX;   // Detection fixed design Matrix
    matrix [n,p_DZ] DZ;   // Detection Random design matrix
    int i_AZ_size;
    int i_DZ_size;
    int i_AZ[i_AZ_size];  // Abundance sigma -> gamma design matrix
    int i_DZ[i_DZ_size];  // Detection sigma -> gamma design matrix
    int N_limit;          // value for summing out the Poisson distribution
  }
  parameters {
    vector <lower=-5, upper=5> [p_AX] beta_A;
    vector <lower=-5, upper=5> [p_AZ] gamma_A;
    vector <lower=-5, upper=5> [p_DX] beta_D;
    vector <lower=-5, upper=5> [p_DZ] gamma_D;
    vector <lower=.1, upper=10> [q_AZ] sigma_A;
    vector <lower=.1, upper=10> [q_DZ] sigma_D;
  }
  transformed parameters{
    vector [n] lambda;
    vector [n] r;
    vector [n] logit_r;
    
    //print("beta_A is:", beta_A);
    //print("beta_D is:", beta_D);

    if( p_AZ > 1 ){
      lambda  <- exp( AX * beta_A + AZ * gamma_A );
      //lambda  <- exp( AX * beta_A );
    }else{
      lambda  <- exp( AX * beta_A );
    }

    if( p_DZ > 1){
      logit_r <- DX * beta_D + DZ * gamma_D;
      //logit_r <- DX * beta_D;
    }else{
      logit_r <- DX * beta_D;
    }
    for(i in 1:n){
      r[i] <- inv_logit(logit_r[i]);
      //lambda[i] <- 1;
      //r[i] <- 0.5;
    }
  }
  model {
    int NextSite;
    real PoissonDensity;
    vector [N_limit+1] logBinomialDensity;

    beta_A ~ normal(0,10);
    beta_D ~ normal(0,10);
    sigma_A ~ uniform(.1,10);
    sigma_D ~ uniform(.1,10);
    for( k in 1:p_AZ ){
      //gamma_A[k] ~ normal(0, 5);
      gamma_A[k] ~ normal(0, sigma_A[ i_AZ[k] ]);
    }
    for( k in 1:p_DZ ){
      //gamma_D[k] ~ normal(0, 5);
      gamma_D[k] ~ normal(0, sigma_D[ i_DZ[k] ]);
    }

    for( N in 0:N_limit ){
      logBinomialDensity[N+1] <- 0;  // Initialize to zero
    }

    for( i in 1:n ){
      if( i == n ){
        NextSite <- Y[i,2] + 1;     // Site number 
      }else{
        NextSite <- Y[i+1,2]; 
      }

      // N = 0      
      if( Y[i,5] >= 1 ){  // if Detections >= 1, N = 0 then cannont be true
        logBinomialDensity[1] <- logBinomialDensity[1] - 10; //+ negative_infinity();
      }else{ // if N = 0, then event Detections = 0 has probability = 1.
        // logBinomialDensity[1] <- logBinomialDensity[1] + 0;
      } 
      for( N in 1:N_limit){
        logBinomialDensity[N+1] <- logBinomialDensity[N+1] + 
                                   binomial_log( Y[i,5], Y[i,4], 1-(1-r[i])^N );
      }

      if( NextSite != Y[i,2] ){  // If we are finished with this site
        PoissonDensity <- 0;
        for( N in 0:N_limit){
          PoissonDensity <- PoissonDensity + 
                exp( logBinomialDensity[N+1] + 
                     poisson_log(N, lambda[i]) );        
        }
        //print("PoissonDensity is ", PoissonDensity);
        increment_log_prob( log(PoissonDensity) );

        for( N in 0:N_limit ){
          logBinomialDensity[N+1] <- 0;  // Reset for the next site
        }
      }

    } // for(i in 1:n)
  } // model

  //generated quantities{
    //vector [n] Psi_fixed;
    //vector [n] Psi_random;
    //vector [n] lambda_fixed;
    //vector [n] lambda_random;

    //lambda_fixed <- exp( AX * beta_A );
    //if( p_AZ > 1 ){
    //  lambda_random  <- exp(AX * beta_A + AZ * gamma_A);
    //}else{
    //  lambda_random <- exp( AX * beta_A);
    //}

    //for(i in 1:n){
      //Psi_fixed[i]  <- 1 - exp( poisson_log(0, lambda_fixed[i]));
      //Psi_random[i] <- 1 - exp( poisson_log(0, lambda_random[i]));
    //}
  //}

