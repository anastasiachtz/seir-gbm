library(rstan)
library(readr)
rstan_options(auto_write = TRUE)           
options(mc.cores = parallel::detectCores())

source("data_to_stan_list.R")
source("initial_values.R")
set.seed(1234)

n_chains=4
n_warmups=2000
n_iter=4500
n_thin=5
n_adapt_delta=0.99
n_max_treedepth=16

nuts_fit <- function(country = c("GR", "PT", "UK", "DE", "SE", "NO"), n_chains, n_warmups, n_iter, n_thin, n_adapt_delta,n_max_treedepth){
  
  # Specify parameters to monitor
  parameters_SEIR = c("beta0", "beta_N", "gamma1", "gamma2_N", 
                      "sigmaBM1", "sigmaBM2","sigmaBM3","sigmaBM4", 
                      "ifr", "phiD",
                      "R_0", "c_tot", "muD",
                      "y_hat", "log_lik", "dev")
  
  fit_by_country <-  stan("SEIR.stan", 
                          data = cov_data(country), 
                          pars = parameters_SEIR , 
                          init = inits(country), 
                          chains = n_chains, 
                          warmup = n_warmups, 
                          iter = n_iter, 
                          thin=n_thin, 
                          control=list(adapt_delta=n_adapt_delta, max_treedepth=n_max_treedepth), seed=4321)
  
  return(fit_by_country)
}

fit_GR <- nuts_fit("GR", n_chains, n_warmups, n_iter, n_thin, n_adapt_delta,n_max_treedepth)
save.image(file='GR_SEIR_VAC.RData')

fit_PT <- nuts_fit("PT", n_chains, n_warmups, n_iter, n_thin, n_adapt_delta,n_max_treedepth)
save.image(file='PT_SEIR_VAC.RData')

fit_UK <- nuts_fit("UK", n_chains, n_warmups, n_iter, n_thin, n_adapt_delta,n_max_treedepth)
save.image(file='UK_SEIR_VAC.RData')

fit_DE <- nuts_fit("DE", n_chains, n_warmups, n_iter, n_thin, n_adapt_delta,n_max_treedepth)
save.image(file='DE_SEIR_VAC.RData')

fit_SE <- nuts_fit("SE", n_chains, n_warmups, n_iter, n_thin, n_adapt_delta,n_max_treedepth)
save.image(file='SE_SEIR_VAC.RData')

fit_NO <- nuts_fit("NO", n_chains, n_warmups, n_iter, n_thin, n_adapt_delta,n_max_treedepth)
save.image(file='NO_SEIR_VAC.RData')
