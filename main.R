library(dplyr)
library(rstan)
library(bayesplot)
library(ggplot2)
library(readr)
rstan_options(auto_write = TRUE)           
options(mc.cores = parallel::detectCores())

##############################################################################
population_UK <- 67886011    #Population numbers: https://www.worldometers.info/world-population/population-by-country/

time_series_covid19_confirmed_global <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")

time_series_covid19_deaths_global <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv")

full_cum_confUK <- subset(time_series_covid19_confirmed_global, `Country/Region`=="United Kingdom")[c(5:length(time_series_covid19_confirmed_global))]
full_cum_deadUK <- subset(time_series_covid19_deaths_global, `Country/Region`=="United Kingdom")[c(5:length(time_series_covid19_deaths_global))]

aggreg_cum_confUK <- as.matrix(full_cum_confUK)
full_cum_confUK <- (as.numeric(colSums(aggreg_cum_confUK)))
aggreg_cum_deadUK <- as.matrix(full_cum_deadUK)
full_cum_deadUK <- (as.numeric(colSums(aggreg_cum_deadUK)))

timeUK <- (seq(as.Date('2020/01/22'),as.Date('2021/09/18'),"days"))
fit_timeUK <- timeUK[38:526]   # 28/2 - 30/06
cum_confUK <- full_cum_confUK[38:526]  
cum_deadUK <- full_cum_deadUK[38:526]
sample_daysUK <- length(fit_timeUK)
sample_time=1:sample_daysUK

new_casesUK <- rep(0,sample_daysUK)
new_casesUK[1] <- cum_confUK[1]
for (t in 2:sample_daysUK){ 
  new_casesUK[t] <- cum_confUK[t]-cum_confUK[t-1]
}
deadUK <- rep(0,sample_daysUK)
deadUK[1] <- cum_deadUK[1]
for (t in 2:sample_daysUK){
  deadUK[t] <- cum_deadUK[t]-cum_deadUK[t-1]
  if (deadUK[t]<0) {
    deadUK[t]=0
  } }

vaccines <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv")
vaccines_UK <- subset(vaccines, `location`=="United Kingdom")
cum_vaccines1 <- vaccines_UK$people_vaccinated
new_vaccines1 <- rep(0,length(cum_vaccines1))
new_vaccines1[1:29] <- c(86465, 0, 0, 0, 0, 0, 0, 588821, 0, 0, 0, 0, 0, 0,  329787, 0, 0, 0, 0, 0, 0, 375357, 0, 0, 0, 0, 0, 0, 906142)
for (t in 30:length(cum_vaccines1)){
  new_vaccines1[t] <- cum_vaccines1[t]-cum_vaccines1[t-1]
}
#length(seq(as.Date('2020/02/28'),as.Date('2020/12/12'),"days")) + 45 
vaccines_t_45days <- c(rep(0,334),new_vaccines1)
vaccines_t_45days <- vaccines_t_45days[1:sample_daysUK]

# Discretize infection to death distribution
infection_death = rep(0,sample_daysUK)
infection_death[1] = pgamma(1.5,shape=6.29,rate =0.26) - pgamma(0,shape=6.29,rate =0.26)
for(i in 2:length(infection_death)) {
  infection_death[i] = pgamma(i+.5,shape=6.29,rate =0.26) - pgamma(i-.5,shape=6.29,rate =0.26)
}

left <- seq(1,sample_daysUK,1)
right <- seq(2,(sample_daysUK+1),1)

# Almost all individuals are susceptible at the start of the epidemic
initial <- c(67884511,450,450,300,300,0,300)

approx_ifr <- c(0.01035,0.007245,0.0095,0.00004,0.000018)

# Modify data into a form suitable for Stan
cov_data_UK = list(n_obs = sample_daysUK,
                n_difeq = 7,
                n_pop = population_UK,
                v = vaccines_t_45days,
                yD = deadUK,
                y_init = initial,
                t0 = 0,
                ts = sample_time,
                I_D = infection_death,
                left_t=left,
                right_t=right,
                time_change_gamma=309,
                rho=0.9,
                sigmaBM_cp1=156,  #1/8
                sigmaBM_cp2=277,  #30/11
                ifr_cp1=142, #18/07
                ifr_cp2=217, #01/10
                ifr_cp3=337, #30/01
                ifr_cp4=460, #01/06
                ifr_mu=approx_ifr) 

# Specify parameters to monitor
parameters = c("beta0", "beta_N", "gamma1", "gamma2_N", "sigmaBM1", "sigmaBM2","sigmaBM3", "ifr", "R_0", "y_hat", "c_tot","muD","phiD", "log_lik", "dev")


set.seed(1234)
# Set initial values:
ini_1 = function(){
  list(eta0= runif(1,0.9,1.1),  
       eta = rep(runif(1,-2.5, 1.1),sample_daysUK),
       gamma1=runif(1,0.95,1.05), 
       gamma2_N=c(runif(1,0.39,0.41), runif(1,0.49,0.51)),
       sigmaBM1=runif(1,0.1,0.8), sigmaBM2=runif(1,0.01,0.1), sigmaBM3=runif(1,10,20),
       ifr=c(runif(1,0.01031,0.01037), runif(1,0.007,0.0075),
             runif(1,0.0093,0.0097), runif(1,0.000038,0.000042),
             runif(1,0.000016,0.00002)),
       reciprocal_phiD=runif(1,0.1,0.5))
}

test = stan("SEIR.stan", data = cov_data_UK, init = ini_1, pars = parameters, chains = 1, iter = 3,control=list(adapt_delta=0.99, max_treedepth=16),seed=4321)

n_chains=4
n_warmups=2000
n_iter=12000
n_thin=10
time.start_nuts1 <- Sys.time()
nuts_fit_1 = stan("SEIR.stan", data = cov_data_UK, pars = parameters, init = ini_1, chains = n_chains, warmup = n_warmups, iter = n_iter, thin=n_thin, control=list(adapt_delta=0.99, max_treedepth=16), seed=4321)
time.end_nuts1 <- Sys.time()
duration_nuts1<- time.end_nuts1 - time.start_nuts1
```