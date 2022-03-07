library(here)
library(rstan)
library(readr)
library(gtools)
library(scales)
rstan_options(auto_write = TRUE)           
options(mc.cores = parallel::detectCores())

source(here::here('regressions',"multi_reg_data_to_stan_list.R"))
source(here::here('regressions',"multi_reg_inits.R"))
set.seed(1234)

n_chains=4
n_warmups=2000
n_iter=6000
n_thin=1
n_adapt_delta=0.99
n_max_treedepth=18

nuts_fit_reg <- function(country = c("GR", "PT", "UK", "DE", "SE", "NO"), n_chains, n_warmups, n_iter, n_thin, n_adapt_delta,n_max_treedepth){
  
  # Specify parameters to monitor
  parameters_regression = c("delta", "L_Omega", "L_sigma", "sigma", 
                            "mu", "L_Sigma", "B0", 
                            "log_lik", "dev", "Sigma", "Omega")
  
  reg_by_country <-  stan("multi_reg.stan", 
                          data = reg_data(country), 
                          pars = parameters_regression , 
                          init = inits_reg(country), 
                          chains = n_chains, 
                          warmup = n_warmups, 
                          iter = n_iter, 
                          thin=n_thin, 
                          control=list(adapt_delta=n_adapt_delta, max_treedepth=n_max_treedepth), seed=4321)
  
  return(fit_reg_by_country)
}

reg_GR <- nuts_fit_reg("GR", n_chains, n_warmups, n_iter, n_thin, n_adapt_delta,n_max_treedepth)
save.image(file='GR_reg.RData')

reg_PT <- nuts_fit_reg("PT", n_chains, n_warmups, n_iter, n_thin, n_adapt_delta,n_max_treedepth)
save.image(file='PT_reg.RData')

reg_UK <- nuts_fit_reg("UK", n_chains, n_warmups, n_iter, n_thin, n_adapt_delta,n_max_treedepth)
save.image(file='UK_reg.RData')

reg_DE <- nuts_fit_reg("DE", n_chains, n_warmups, n_iter, n_thin, n_adapt_delta,n_max_treedepth)
save.image(file='DE_reg.RData')

reg_SE <- nuts_fit_reg("SE", n_chains, n_warmups, n_iter, n_thin, n_adapt_delta,n_max_treedepth)
save.image(file='SE_reg.RData')

reg_NO <- nuts_fit_reg("NO", n_chains, n_warmups, n_iter, n_thin, n_adapt_delta,n_max_treedepth)
save.image(file='NO_reg.RData')
