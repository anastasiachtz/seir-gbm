library(deSolve)
library(rstan)
rstan_options(auto_write = TRUE)           
options(mc.cores = parallel::detectCores())

set.seed(1234)

#Simulate data

# Define time-varying sigma function
time_varying_sigma <- function(time) {
  if (time < 90) {
    return(0.1)
  } else {
    return(0.05)
  }
}
# Function to generate geometric Brownian motion for beta
generate_gbm <- function(initial_beta, mu, sigma, times) {
  dt <- diff(times)[1]
  n <- length(times)
  
  beta <- numeric(n)
  beta[1] <- initial_beta
  
  for (t in 2:n) {
    sigma <- time_varying_sigma(times[t])
    dW <- rnorm(1, mean = 0, sd = sqrt(dt))
    beta[t] <- beta[t-1] * exp((mu - 0.5 * sigma^2) * dt + sigma * dW)
  }
  
  return(beta)
}

# SEIR model function with time-varying beta
seir_model_gbm <- function(time, state, parameters, beta_values) {
  with(as.list(c(state, parameters)), {

    beta <- beta_values[which.min(abs(time - times))]
    gamma1 <- parameters["gamma1"]
    gamma2 <- parameters["gamma2"]
    
    population <- 100000
  
    dS <- - beta * S * (I1 + I2)/population
    dE1 <- (beta * S * (I1 + I2)/population) - gamma1 * E1
    dE2 <- gamma1 * (E1 - E2)
    dI1 <- gamma1 * E2 - gamma2 * I1
    dI2 <- gamma2 * (I1 - I2)
    dR <- gamma2 * I2
    dC <- gamma1 * E2  # Dummy compartment for new cases
    
    # Return the derivatives
    return(list(c(dS, dE1, dE2, dI1, dI2, dR, dC)))
  })
}

initial_state <- c(S = 99660, E1 = 100, E2 = 100, 
                   I1 = 70, I2 = 70, R = 0, C = 0) 

# Parameters
parameters <- c(gamma1 = 2/3, gamma2 = 2/5)

# Time vector
times <- seq(0, 180, by = 1)

set.seed(123)
# Generate the beta values using GBM
initial_beta <- 0.3 # initial value for beta
mu <- 0         # drift coefficient
#sigma <- 0.5   # volatility coefficient

beta_values <- generate_gbm(initial_beta, mu, sigma, times)
beta_df <- data.frame(time = times, beta = beta_values)

# Solve the differential equations
out <- ode(y = initial_state, times = times, func = seir_model_gbm, parms = parameters, beta_values = beta_values)


out_df <- as.data.frame(out)

# Calculate daily new cases
out_df$new_cases <- c(NA, diff(out_df$C))

n_obs <- length(times)-1  # Number of observations
ifr <- 0.1                # Infection Fatality Rate (IFR)
infection_death = rep(0,n_obs) # Discretized infection to death distribution
infection_death[1] = pgamma(1.5,shape=6.29,rate =0.26) - pgamma(0,shape=6.29,rate =0.26)
for(i in 2:n_obs) {
  infection_death[i] = pgamma(i+.5,shape=6.29,rate =0.26) - pgamma(i-.5,shape=6.29,rate =0.26)
}

infection_death_rev <- rev(infection_death)  # Example I_D_rev (reversed c_tot)

# Calculate daily new deaths
deaths <- numeric(n_obs)
deaths[1] <- 0
for (i in 2:n_obs) {
  deaths[i] <- ifr * sum(head(out_df$new_cases[-1], i - 1) * tail(infection_death_rev, i - 1))
}
deaths <- c(NA, deaths)
# Apply Gaussian kernel smoothing to deaths, ignoring the first NA value
gaussian_smooth <- function(x, bandwidth = 5) {
  ksmooth(1:length(x), x, "normal", bandwidth = bandwidth)$y
}

smoothed_deaths <- gaussian_smooth(deaths[-1], bandwidth = 7)  # 7-day bandwidth
smoothed_deaths <- c(NA, smoothed_deaths)
smoothed_deaths <- round(smoothed_deaths) # Convert deaths to integers
out_df$deaths <- smoothed_deaths

############################################################################################
# Plot the results
# library(ggplot2)
# 
# # Plot the beta values over time
# ggplot(beta_df, aes(x = time, y = beta)) +
#   geom_line(color = "blue") +
#   labs(title = "Time-varying Beta (Geometric Brownian Motion)", x = "Time", y = "Beta") +
#   theme_minimal()
# 
# ggplot(out_df, aes(x = time)) +
#   geom_line(aes(y = S, color = "Susceptible")) +
#   geom_line(aes(y = E1, color = "Exposed1")) +
#   geom_line(aes(y = E2, color = "Exposed2")) +
#   geom_line(aes(y = I1, color = "Infectious1")) +
#   geom_line(aes(y = I2, color = "Infectious2")) +
#   geom_line(aes(y = R, color = "Recovered")) +
#   labs(y = "Population", color = "Compartment") +
#   theme_minimal()
# 
# # Plot the new cases over time
# ggplot(out_df, aes(x = time, y = new_cases)) +
#   geom_line(color = "red") +
#   labs(title = "New Cases Over Time", x = "Time", y = "New Cases") +
#   theme_minimal()
# 
# # Plot the deaths over time
# ggplot(out_df, aes(x = time, y = deaths)) +
#   geom_line(color = "black") +
#   labs(title = "Deaths Over Time", x = "Time", y = "Deaths") +
#   theme_minimal()

######################################################################
left <- seq(1,length(times)-1,1)
right <- seq(2,length(times),1)

sim_data = list(n_obs = length(times)-1,
                           n_difeq = 7,
                           n_pop = 100000,
                           yD = out_df$deaths[-1],
                           y_init = initial_state,
                           t0 = 0,
                           ts = times[-1],
                           I_D = infection_death,
                           left_t=left,
                           right_t=right,
                           sigmaBM_cp=90,
                           ifr_mu=0.1)
########################################################################
n_chains=3
n_warmups=500
n_iter=2500
n_thin=2
n_adapt_delta=0.85
n_max_treedepth=12

# Specify parameters to monitor
parameters_SEIR = c("beta0", "beta_N", "gamma1", "gamma2", 
                      "sigmaBM1", "sigmaBM2", 
                      "ifr", "phiD",
                      "R_0", "c_tot", "muD",
                      "y_hat", "log_lik", "dev")


inits <- function(){
          list(eta0= runif(1,-1.5,-0.5),
           eta = c(rep(runif(1,-1.5,-0.5),length(times)-1)),
           gamma1=runif(1,0.8,1.2), 
           gamma2=runif(1,0.3,0.5),
           sigmaBM1=runif(1,0.01,0.2), sigmaBM2=runif(1,0.01,0.2),
           ifr=runif(1,0.098,0.11),
           reciprocal_phiD=runif(1,0.005,0.009))
}
time.start_nuts <- Sys.time()  
fit <-  stan("SEIR.stan", 
        data = sim_data, 
        pars = parameters_SEIR , 
        init = inits, 
        chains = n_chains, 
        warmup = n_warmups, 
        iter = n_iter, 
        thin=n_thin, 
        control=list(adapt_delta=n_adapt_delta, max_treedepth=n_max_treedepth), 
        seed=4321)
time.end_nuts <- Sys.time()
duration_nuts<- time.end_nuts - time.start_nuts
save.image(file='fit1.RData')
