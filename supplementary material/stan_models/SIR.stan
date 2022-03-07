functions {
  real[] SIR(real time,          // time
                real[] y,           // system state {susceptible,infected,recovered}
                real[] theta,       // parameters
                real[] x_r,
                int[] x_i) {
    // system states           
    real S = y[1];
    real I1 = y[2];
    real I2 = y[3];
    real R = y[4];
    
    int n_pop = x_i[1];    
    int n_obs = x_i[2];
    int v[n_obs] = x_i[3:n_obs+2]; 
    real left_t[n_obs] = x_r[1:n_obs];
    real right_t[n_obs] = x_r[n_obs+1:2*n_obs];
    real rho = x_r[2*n_obs+1];
    real time_change_gamma = x_r[2*n_obs+2];
    
    real dy_dt[5];
    
    // estimated parameters
    real beta0 = theta[1];
    real beta_N[n_obs] = theta[2:(n_obs+1)];
    real gamma_N[2] = theta[(n_obs+2):(n_obs+3)];
    
    real beta;
    real gamma;
    real vaccinated;
    
    for (i in 1:n_obs) {
     if(time >= 0 && time < 1){
        beta = beta0;
        vaccinated = rho * 0;
     }
      else if(time >= left_t[i] && time < right_t[i]){
      beta = beta_N[i];
      vaccinated = rho * v[i];
      }
    }
    
    if(time < time_change_gamma){
       gamma = gamma_N[1];
    }
    else if(time >= time_change_gamma){
       gamma = gamma_N[2];
    }
    
    // SIIR in Numbers (NOT fractions of the population)
    //S
    dy_dt[1] = - (beta * S * (I1+I2) / n_pop) - vaccinated;
    //I1
    dy_dt[2] = (beta * S * (I1+I2) / n_pop) - (gamma * I1);
    //I2
    dy_dt[3] = (gamma * I1) - (gamma * I2);
    //R
    dy_dt[4] = (gamma * I2)  + vaccinated;
    //C
    dy_dt[5] = beta * S * (I1+I2) / n_pop;  // add a dummy compartment to record the cumulative incidence
    
    
    return dy_dt;
  }
  
  real[ , ] integrate_ode_trapezoidal(real[] y_initial, real initial_time, real[] times, real[] theta, real[] x_r, int[] x_i) {
    real h;
    vector[size(y_initial)] dy_dt_initial_time;
    vector[size(y_initial)] dy_dt_t;
    vector[size(y_initial)] k;
    real y_approx[size(times),size(y_initial)];
    
    h = times[1] - initial_time;
    dy_dt_initial_time = to_vector(SIR(initial_time, y_initial, theta, x_r, x_i));
    k = h*dy_dt_initial_time;
    y_approx[1,] = to_array_1d(to_vector(y_initial) + h*(dy_dt_initial_time + 
                                                           to_vector(SIR(times[1], to_array_1d(to_vector(y_initial)+k), theta, x_r, x_i)))/2);
    
    for (t in 1:size(times)-1) {
      h = (times[t+1] - times[t]);
      dy_dt_t = to_vector(SIR(times[t], y_approx[t], theta, x_r, x_i));
      k = h*dy_dt_t;
      y_approx[t+1,] = to_array_1d(to_vector(y_approx[t,]) + h*(dy_dt_t + 
                                                                  to_vector(SIR(times[t+1], to_array_1d(to_vector(y_approx[t,])+k), theta, x_r, x_i)))/2);
    }
    
    return y_approx;
  }
}

data {
  int<lower = 1> n_obs;         // number of days observed
  int<lower = 1> n_difeq;       // number of differential equations
  int<lower = 1> n_pop;         // population 
  int<lower = 0> v[n_obs];      // data, daily number of vaccinations          
  int yD[n_obs];                // data, daily number of deaths
  int<lower = 1> sigmaBM_cp1;   // change points on BM volatility
  int<lower = 1> sigmaBM_cp2;   
  int<lower = 1> sigmaBM_cp3;
  int<lower = 1> ifr_cp1;       // change points on ifr
  int<lower = 1> ifr_cp2;
  int<lower = 1> ifr_cp3;
  int<lower = 1> ifr_cp4;
  real y_init[n_difeq];         // initial conditions for S,E,I,R,C
  real t0;                      // initial time point (zero)
  real ts[n_obs];               // time points observed
  real<lower=0> left_t[n_obs];  // left time limit
  real<lower=0> right_t[n_obs]; // right time limit
  real<lower=0> time_change_gamma;  // right time limit for change in gamma
  real I_D[n_obs];              // discretized infection to death distribution
  real<lower = 0,upper=1> rho;  // efficacy
  real<lower=0> ifr_mu[5];      // mean ifr
}

transformed data {
  real x_r[2*n_obs+2];
  int x_i[n_obs+2];
  real I_D_rev[n_obs];    // reversed discretized infection to death distribution
  
  x_i[1] = n_pop;
  x_i[2] = n_obs;
  x_i[3:n_obs+2] = v;
  x_r[1:n_obs] = left_t;
  x_r[n_obs+1:2*n_obs] = right_t;
  x_r[2*n_obs+1] = rho;
  x_r[2*n_obs+2] = time_change_gamma;
  
  for(i in 1:n_obs){
    I_D_rev[i] = I_D[n_obs-i+1];
  }
}

parameters {
  real eta0;                         // initial BM
  real eta[n_obs];                   // BM                     
  real<lower = 0> gamma_N[2];       // recovery rate
  real<lower = 0> sigmaBM1;          // standard deviation of BM 1st wave
  real<lower = 0> sigmaBM2;          // standard deviation of BM 2nd wave
  real<lower = 0> sigmaBM3;          // standard deviation of BM 3rd wave
  real<lower = 0> sigmaBM4;          // standard deviation of BM 4th wave
  real<lower = 0> ifr[5];            // probability of death given infection
  real<lower = 0> reciprocal_phiD;   // overdispersion of NB on deaths
}

transformed parameters{
  real<lower = 0> beta0;                // initial transmission rate
  real<lower = 0> beta_N[n_obs];        // transmission rate, beta=exp(z)
  real theta[n_obs+3];                  // ode parameters vector
  real y_hat[n_obs, n_difeq];           // solution from the ODE solver 
  real<lower = 0> c_tot[n_obs];         // total new cases
  real<lower = 0> muD[n_obs];           // mean of NB on deaths
  real<lower = 0> phiD;                 // 1/reciprocal_phiD
  
  beta0 = exp(eta0);
  for (i in 1:n_obs){
    beta_N[i] = exp(eta[i]);
  }
  
  theta[1] = beta0;
  theta[2:(n_obs+1)] = beta_N;
  theta[(n_obs+2):(n_obs+3)] = gamma_N;
  
  //y_hat = integrate_ode_rk45(SIR, y_init, t0, ts, theta, x_r, x_i,1e-10, 1e-10, 10000);
  //y_hat = integrate_ode_bdf(SIR, y_init, t0, ts, theta, x_r, x_i,1e-10, 1e-10, 10000);
  y_hat = integrate_ode_trapezoidal(y_init, t0, ts, theta, x_r, x_i);
  
  for (i in 1:n_obs){
    c_tot[i]=y_hat[i,4]- (i==1 ? 0 : (y_hat[i,4]>y_hat[i-1,4] ? y_hat[i-1,4] : 0) );
  }
  
  muD[1] = 1e-010; 
  for (i in 2:n_obs){
    if (i<ifr_cp1)
      muD[i] =  ifr[1] * dot_product(head(c_tot,i-1),tail(I_D_rev, i-1));
    else if (ifr_cp1<=i<ifr_cp2)
      muD[i] =  ifr[2] * dot_product(head(c_tot,i-1),tail(I_D_rev, i-1));
    else if (ifr_cp2<=i<ifr_cp3)
      muD[i] =  ifr[3] * dot_product(head(c_tot,i-1),tail(I_D_rev, i-1));
    else if (ifr_cp3<=i<ifr_cp4)
      muD[i] =  ifr[4] * dot_product(head(c_tot,i-1),tail(I_D_rev, i-1));
    else
      muD[i] =  ifr[5] * dot_product(head(c_tot,i-1),tail(I_D_rev, i-1));
  }
  
  phiD = 1. / reciprocal_phiD;
}

model {
  //priors
  eta0 ~ normal(0, 1);
  gamma_N[1] ~ gamma(1111,3333);    // Assume mean infectious period = 6 days, mean=2/6
  gamma_N[2] ~ gamma(1600,4000);    // Assume mean infectious period = 5 days, mean=2/5
  sigmaBM1 ~ cauchy(0,5);
  sigmaBM2 ~ cauchy(0,5);
  sigmaBM3 ~ cauchy(0,5);
  sigmaBM4 ~ cauchy(0,5);
  for (i in 1:5){
    ifr[i] ~ beta_proportion(ifr_mu[i],10000000000);
  }
  reciprocal_phiD ~ cauchy(0,5);
  
  eta[1] ~ normal(eta0, sigmaBM1);    // centered parameterization
  for (i in 2:n_obs){
    if (i<sigmaBM_cp1)
      eta[i] ~ normal(eta[i-1], sigmaBM1);
    else if (sigmaBM_cp1<=i<sigmaBM_cp2)
      eta[i] ~ normal(eta[i-1], sigmaBM2);
    else if (sigmaBM_cp2<=i<sigmaBM_cp3)
      eta[i] ~ normal(eta[i-1], sigmaBM3);
    else
      eta[i] ~ normal(eta[i-1], sigmaBM4);
  }
  
  yD ~ neg_binomial_2(muD,phiD);
}

generated quantities {
  vector[n_obs] R_0;           // Basic reproduction number
  vector[n_obs] log_lik;       // log-likelihood for loo package
  real dev;                    // deviance
  
  for (i in 1:n_obs){
  if (i<time_change_gamma)
  R_0[i] = (2*beta_N[i])/gamma_N[1];
  else
  R_0[i] = (2*beta_N[i])/gamma_N[2];
  }
  
  for (i in 1:n_obs) {
    log_lik[i] = neg_binomial_2_lpmf(yD[i]| muD[i],phiD);
  }
  
  dev = (-2) * sum(log_lik);
}