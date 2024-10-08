
# A modelling framework for the analysis of the transmission of SARS-CoV2

## Citation

A. Chatzilena, N. Demiris, and K. Kalogeropoulos. A modelling framework for the analysis of the transmission of
SARS-CoV2. Statistics in Medicine. 2024; 1-17. [https://doi.org/10.1002/sim.10195](https://doi.org/10.1002/sim.10195)
### Background

Action plans against the current SARS-CoV-2 pandemic have been implemented around the globe in order to reduce transmission. The reproduction number has been found to respond to public health interventions changing throughout the pandemic waves. However, the actual global burden of SARS-CoV-2 remains unknown due to severe under-ascertainment of cases. The use of reported  deaths has been pointed out as a more reliable source of information, likely less prone to under-reporting. Given that daily deaths occur from past infections weighted by their probability of death, one may infer the true number of infections accounting for their age distribution, using the data on reported deaths.

### Methods 

This work uses a model-based approach to estimate the transmissibility of SARS-CoV-2 and the effect of the adopted control measures across 6 European countries. We extend the [Dureau, Kalogeropoulos and Baguelin (2013)](https://doi.org/10.1093/biostatistics/kxs052) model by introducing an additional hierarchical level to capture the unknown true number of cases through deaths and we fit to the unknown true number of cases an SEIR (susceptible-exposed-infected-recovered) compartmental model driven by a stochastic time-varying transmission rate that captures the effect of both the control measures and the behavioural changes. Our model may also be viewed as an extension of the work in [Flaxman et al. (2020)](https://doi.org/10.1038/s41586-020-2405-7) where instead of employing a transmission model with instantaneous intervention effects on the reproduction number, the estimated number of cases is indirectly inferred and generated by a diffusion-driven stochastic SEIR process.

We implement our suggested model in a Bayesian framework, using Hamiltonian Monte Carlo employing the Stan software, by fitting our model to daily reported deaths for Greece, Portugal, United Kingdom, Germany, Sweden and Norway. We then examine how we can combine our estimates of the total number of daily cases with data on daily laboratory-confirmed cases and estimate the daily reporting ratio. Finally, through a multivariate regression analysis, we disentangle the effects of preventive measures and testing policies on the estimated total cases, the time-varying transmission rate and the reporting rate, using only publicly available data. Our proposed model is then fit to data from country pairs with similar population demographics and health and social welfare infrastructures, to gain insights on the COVID-19 pandemic.

### Results

We estimate that during the time course of the pandemic there have been substantially more infections than those detected by health care systems. Especially during the peak of each pandemic wave, the actual number of infections is significantly larger. The estimated cumulative cases can offer a measure of the actual burden of the pandemic. We show that the estimated changes in the reproduction number are consistent with the expected variation in SARS-CoV-2 transmission over time, as a result of the implemented control strategies. We estimate that all countries except Sweden, having introduced several non-pharmaceutical interventions, were able to drop Rt below 1 well short of their nationwide lockdowns. The effects of sequentially introduced interventions in a small period of time, are highly interdependent making it hard to disentangle their individual contribution. Therefore, one may not offer robust conclusions concerning optimal strategies in the absence of additional data sources.

### Conclusions

A distinct characteristic of our modelling approach is the absence of strong structural assumptions for the temporal evolution of the transmission rate, and subsequently for the case reproduction number, departing from the piecewise constant assumption of  [Flaxman et al. (2020)](https://doi.org/10.1038/s41586-020-2405-7). Changes in Rt are only driven by variations in the observed data on deaths for each country. We consider that control measures, different public responses on these measures based on cultural characteristics, adaptive human behaviour during a pandemic and any other time-varying factor is reflected in the trends in the numbers of deaths resulting from the respective infections. Therefore, our framework can be adapted to other countries in a straightforward manner.

Several limitations need to be considered when using death counts as the main source of information. Early in the pandemic, in the absence of European and international standards, some deaths due to COVID-19 may not be recorded leading to underestimation of infections. Therefore our initial estimates must be viewed with caution. Also, reporting procedures differ between countries both in terms of the timing of the report as well as the definition of COVID-19 related death. In addition, data on deaths depend on past infections and are not suitable for real-time analysis without further assumptions. However, in the absence of large seroprevalence studies in many countries, death counts offer a credible option for evaluating the actual burden of the pandemic in terms of people infected.

# Code

The following files/folders are available.

[main.R](https://github.com/anastasiachtz/seir-gbm/blob/main/main.R) runs the analysis for each country using Stan’s No-U-Turn sampler variant of Hamiltonian Monte Carlo and saves the outputs.

[SEIR.stan](https://github.com/anastasiachtz/seir-gbm/blob/main/SEIR.stan) contains the [Stan](https://mc-stan.org/) model statement where daily new infections are generated by a stochastic transmission SEIR compartmental model and linked to the reported deaths.

[data_to_stan_list.R](https://github.com/anastasiachtz/seir-gbm/blob/main/data_to_stan_list.R) contains a function which modifies the data into a form suitable for Stan for each country.

[initial_values.R](https://github.com/anastasiachtz/seir-gbm/blob/main/initial_values.R) contains a function which provides initial values for the model parameters for each country.

[data](https://github.com/anastasiachtz/seir-gbm/tree/main/data):<br>
Contains input data that are not directly available from external links.

[regressions](https://github.com/anastasiachtz/seir-gbm/tree/main/regressions):<br>
[main_regressions.R](https://github.com/anastasiachtz/seir-gbm/blob/main/regressions/main_regressions.R) runs the the multivariate regression analysis for each country using Stan’s No-U-Turn sampler variant of Hamiltonian Monte Carlo, [multi_reg.stan](https://github.com/anastasiachtz/seir-gbm/blob/main/regressions/multi_reg.stan) contains the [Stan](https://mc-stan.org/) multivariate regression model, [multi_reg_data_to_stan_list.R](https://github.com/anastasiachtz/seir-gbm/blob/main/regressions/multi_reg_data_to_stan_list.R) modifies the model's posterior estimates and data for the regression, into a form suitable for Stan, 
[multi_reg_inits.R](https://github.com/anastasiachtz/seir-gbm/blob/main/regressions/multi_reg_inits.R) provides initial values for the regression parameters.

[supplementary material](https://github.com/anastasiachtz/seir-gbm/tree/main/supplementary%20material):<br>
[Extended_Supplement_A_modelling_framework_for_the_analysis_of_the_transmission_of_SARS_CoV2.pdf](https://github.com/anastasiachtz/seir-gbm/blob/main/supplementary%20material/Extended_Supplement_A_modelling_framework_for_the_analysis_of_the_transmission_of_SARS_CoV2.pdf) includes implementation details, supplemental tables, figures and diagnostics for the proposed SEIR model, as well as indicative figures for the SEIR model without vaccinations and the analogous SIR model with vaccinations. Analytical results on the posterior estimates of the multivariate regression coefficients and detailed data references.<br>
[stan_models](https://github.com/anastasiachtz/seir-gbm/tree/main/supplementary%20material/stan_models) contains the [Stan](https://mc-stan.org/) model statements for the auxiliary models used [SIR.stan](https://github.com/anastasiachtz/seir-gbm/blob/main/supplementary%20material/stan_models/SIR.stan), [SEIR_noVac.stan](https://github.com/anastasiachtz/seir-gbm/blob/main/supplementary%20material/stan_models/SEIR_noVac.stan).
[simulation_gbm_D.R](https://github.com/anastasiachtz/seir-gbm/blob/main/supplementary%20material/simulation_gbm_D.R) runs a simulation experiment.

