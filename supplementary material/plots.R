library(bayesplot)
library(ggplot2)

# Diagnostics

nuts_fit_summary <- summary(fit, pars = c("lp__","beta0", "beta_N[20]", "beta_N[60]", "beta_N[100]", "beta_N[120]", "beta_N[180]", "gamma1", "gamma2", "sigmaBM1", "sigmaBM2", "ifr", "phiD"))$summary
options(scipen = 999)
print(nuts_fit_1_summary,scientific=FALSE,digits=2)

mod_diagnostics <-rstan::get_sampler_params(fit)

# Check for divergent transitions
check_divergences(fit)
check_treedepth(fit)

posterior <- as.array(fit)
color_scheme_set("viridis")
# Markov chain traceplots
mcmc_trace(posterior, pars=c("lp__","beta0", "beta_N[20]", "beta_N[60]", "beta_N[100]", "beta_N[120]", "beta_N[180]", "gamma1", "gamma2", "sigmaBM1", "sigmaBM2", "ifr", "phiD"))

# Univariate and bivariate marginal posterior distributions
pairs(fit, pars = c("beta0", "beta_N[20]", "beta_N[60]", "beta_N[100]", "beta_N[120]", "beta_N[180]", "gamma1", "gamma2", "sigmaBM1", "sigmaBM2", "ifr", "phiD"), labels = c("beta0", "beta_N[20]", "beta_N[60]", "beta_N[100]", "beta_N[120]", "beta_N[180]", "gamma1", "gamma2", "sigmaBM1", "sigmaBM2", "ifr", "phiD"), 
      cex.labels=1.5, font.labels=9, condition = "accept_stat__")

# Kernel density estimates of each Markov chain separately, overlaid
mcmc_dens_overlay(posterior, pars=c("lp__","beta0", "beta_N[20]", "beta_N[60]", "beta_N[100]", "beta_N[120]", "beta_N[180]", "gamma1", "gamma2", "sigmaBM1", "sigmaBM2", "ifr", "phiD"))

###########################################################################################################################
posts <-  rstan::extract(fit)

post_beta <- posts$beta_N
median_beta = apply(post_beta, 2, median)
low_beta = apply(post_beta, 2, quantile, probs=c(0.025))
high_beta = apply(post_beta, 2, quantile, probs=c(0.975))

fit_cases <- posts$c_tot
median_fit_cases = apply(fit_cases, 2, median)
low_fit_cases = apply(fit_cases, 2, quantile, probs=c(0.025))
high_fit_cases = apply(fit_cases, 2, quantile, probs=c(0.975))

fit_dead <- posts$muD
median_fit_dead = apply(fit_dead, 2, median)
low_fit_dead = apply(fit_dead, 2, quantile, probs=c(0.025))
high_fit_dead = apply(fit_dead, 2, quantile, probs=c(0.975))

# Plot fit deaths
fitD <- ggplot() +
  geom_ribbon(data = data.frame(median_fit_dead, low_fit_dead, high_fit_dead, times[-1]),aes(x=times[-1], ymin = low_fit_dead, ymax = high_fit_dead, fill = "95% CI Fitted deaths"), alpha = 0.6, show.legend = FALSE) +
  geom_line(data = data.frame(median_fit_dead, low_fit_dead, high_fit_dead, times[-1]), aes(x=times[-1], y=median_fit_dead, color = "Median Fitted deaths"), size = 3, show.legend = FALSE) +
  geom_point(data = data.frame(out_df$deaths[-1], times[-1]), aes(x=times[-1], y=out_df$deaths[-1], color="Data"), shape = 19, size = 8, show.legend = FALSE) +
  scale_fill_manual(values = c("95% CI Fitted deaths"="pink4")) +
  scale_colour_manual(name='', values=c('Data'='black','Median Fitted deaths'='salmon4'))+
  labs(x = "Time", y = "Number of New Deaths")+
    guides(colour=guide_legend(override.aes=list(shape=c(16,NA), linetype=c(0,1))))+
   scale_x_continuous(breaks = seq(0, max(times[-1]), by = 20))+
   theme_bw(base_size=65)+
   theme(axis.text.x=element_text(size=55, margin = margin(t = 60)))+
   theme(axis.text.y=element_text(size=55))+
   #theme(legend.position=c(.1, 0.75))+
   #theme(legend.title=element_blank())
   theme(legend.position = "none")
 ggsave("death_plot.png", plot = fitD, width = 3000, height = 2200, units = "px", dpi = 72, limitsize = FALSE) 
 
 
 # Plot fit cases
 fitC <- ggplot() +
   geom_ribbon(data = data.frame(median_fit_cases, low_fit_cases, high_fit_cases, times[-1]),aes(x=times[-1], ymin = low_fit_cases, ymax = high_fit_cases, fill = "95% CI"), alpha = 0.6, show.legend = FALSE) +
   geom_line(data = data.frame(median_fit_cases, low_fit_cases, high_fit_cases, times[-1]), aes(x=times[-1], y=median_fit_cases, color = "Median"), size = 3, show.legend = FALSE) +
   scale_fill_manual(name='', values = c("95% CI"="skyblue3")) +
   scale_colour_manual(name='', values=c('Median'='navyblue','Simulated Cases'='black'))+
   labs(x = "Time", y = "Number of New Cases") +
   geom_point(data = data.frame(out_df$new_cases[-1],times[-1]), aes(x=times[-1], y=out_df$new_cases[-1], color="Simulated Cases"), shape = 19, size = 8, show.legend = FALSE) +
   guides(colour=guide_legend(override.aes=list(shape=c(NA,16), linetype=c(1,0))))+
   scale_x_continuous(breaks = seq(0, max(times[-1]), by = 20))+
   theme_bw(base_size=65)+
   theme(axis.text.x=element_text(size=55, margin = margin(t = 60)))+
   theme(axis.text.y=element_text(size=55))+
   #theme(legend.position=c(.1, 0.75))+
   #theme(legend.title=element_blank())
   theme(legend.position = "none")
 ggsave("case_plot.png", plot = fitC, width = 3000, height = 2200, units = "px", dpi = 72, limitsize = FALSE) 
  
# Plot transmission rate 
  bt <- ggplot() +
    geom_ribbon(data = data.frame(median_beta, low_beta, high_beta, times[-1]),aes(x=times[-1], ymin = low_beta, ymax = high_beta, fill = "95% CI"), alpha = 0.6, show.legend = FALSE)+
    geom_line(data = data.frame(median_beta, low_beta, high_beta, times[-1]), aes(x=times[-1], y=median_beta, color = "Median"), size = 3, show.legend = FALSE) +
    geom_point(data = data.frame(beta_df$beta[-1], times[-1]), aes(x=times[-1], y=beta_df$beta[-1], color="Simulated betas"), shape = 19, size = 8, show.legend = FALSE) +
    scale_colour_manual(name='', values=c('Simulated betas'='black', 'Median'='blue'))+
    scale_fill_manual(values = c("95% CI"="steelblue1")) +
    guides(colour = guide_legend(override.aes = list(shape=c(20),  linetype=c(0))))+
    labs(x = "Time", y = expression(beta[t]))+
    scale_x_continuous(breaks = seq(0, max(times[-1]), by = 20))+
    theme_bw(base_size=65)+
    theme(axis.text.x=element_text(size=55, margin = margin(t = 60)))+
    theme(axis.text.y=element_text(size=55))+
    #theme(legend.position=c(.1, 0.75))+
    #theme(legend.title=element_blank())
    theme(legend.position = "none")
  
  ggsave("beta_plot.png", plot = bt, width = 4000, height = 2200, units = "px", dpi = 72, limitsize = FALSE) 