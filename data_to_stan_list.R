cov_data <- function(country = c("GR", "PT", "UK", "DE", "SE", "NO")){
  
  if (country == "GR"){
    population_GR <- 10423054    #Population numbers: https://www.worldometers.info/world-population/population-by-country/
    time_series_covid19_confirmed_global <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
    time_series_covid19_deaths_global <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv")
    full_cum_confGR <- subset(time_series_covid19_confirmed_global, `Country/Region`=="Greece")[c(5:length(time_series_covid19_confirmed_global))]
    full_cum_deadGR <- subset(time_series_covid19_deaths_global, `Country/Region`=="Greece")[c(5:length(time_series_covid19_deaths_global))]
    timeGR <- as.vector(names(full_cum_confGR))
    cum_confGR <- as.numeric(full_cum_confGR)[(match('3/7/20', timeGR)):(match('6/30/21', timeGR))]
    cum_deadGR <- as.numeric(full_cum_deadGR)[(match('3/7/20', timeGR)):(match('6/30/21', timeGR))]
    fit_timeGR <- timeGR[(match('3/7/20', timeGR)):(match('6/30/21', timeGR))]   # 07/03 - 30/06
    sample_daysGR <- length(fit_timeGR)
    sample_time=1:sample_daysGR
    new_casesGR <- rep(0,sample_daysGR)
    new_casesGR[1] <- cum_confGR[1]
    for (t in 2:sample_daysGR){ 
      new_casesGR[t] <- cum_confGR[t]-cum_confGR[t-1]
    }
    deadGR <- rep(0,sample_daysGR)
    deadGR[1] <- cum_deadGR[1]
    for (t in 2:sample_daysGR){
      deadGR[t] <- cum_deadGR[t]-cum_deadGR[t-1]
      if (deadGR[t]<0) {
        deadGR[t]=0
      } 
    }
    
    Greece_data_V <- read.csv('./data/Greece_vaccinations.csv', head = TRUE, sep=";")
    cum_vaccines1 <- Greece_data_V$X1_dose
    new_vaccines1 <- rep(0,length(cum_vaccines1))
    new_vaccines1[1] <- cum_vaccines1[1]
    for (t in 2:length(cum_vaccines1)){
      new_vaccines1[t] <- cum_vaccines1[t]-cum_vaccines1[t-1]
    }
    vaccines_t_45days <- c(rep(0,(length(seq(as.Date('2020/03/07'),as.Date('2020/12/26'),"days")) + 45)),new_vaccines1)
    vaccines_t_45days_GR <- vaccines_t_45days[1:sample_daysGR]
    
    # Discretize infection to death distribution
    infection_death = rep(0,sample_daysGR)
    infection_death[1] = pgamma(1.5,shape=6.29,rate =0.26) - pgamma(0,shape=6.29,rate =0.26)
    for(i in 2:sample_daysGR) {
      infection_death[i] = pgamma(i+.5,shape=6.29,rate =0.26) - pgamma(i-.5,shape=6.29,rate =0.26)
    }
    
    # Almost all individuals are susceptible at the start of the epidemic
    initial_GR <- c(10422884,50,50,35,35,0,65)
    
    left <- seq(1,sample_daysGR,1)
    right <- seq(2,(sample_daysGR+1),1)
    approx_ifr_GR <- c(0.01142,0.007994,0.01142,0.0115,0.0002)
    
    # Modify data into a form suitable for Stan
    cov_data_by_country = list(n_obs = sample_daysGR,
                    n_difeq = 7,
                    n_pop = population_GR,
                    v = vaccines_t_45days_GR,
                    rho=0.9,
                    yD = deadGR,
                    y_init = initial_GR,
                    t0 = 0,
                    ts = sample_time,
                    I_D = infection_death,
                    left_t=left,
                    right_t=right,   
                    time_change_gamma=match('1/1/21', fit_timeGR),
                    sigmaBM_cp1=match('8/1/20', fit_timeGR),
                    sigmaBM_cp2=match('1/10/21', fit_timeGR),
                    ifr_cp1=match('7/28/20', fit_timeGR),
                    ifr_cp2=match('11/1/20', fit_timeGR),
                    ifr_cp3=match('2/9/21', fit_timeGR),
                    ifr_cp4=match('5/9/21', fit_timeGR),
                    ifr_mu=approx_ifr_GR) 
  }
  
  else if (country == "PT"){
    population_PT <- 10196707
    portugal_data <- read_csv("https://raw.githubusercontent.com/dssg-pt/covid19pt-data/master/data.csv")
    full_cum_confPT <- portugal_data$confirmados
    full_cum_deadPT <- portugal_data$obitos
    timePT <- portugal_data$data
    cum_confPT <- as.numeric(full_cum_confPT)[(match('07-03-2020', timePT)):(match('30-06-2021', timePT))]
    cum_deadPT <- as.numeric(full_cum_deadPT)[(match('07-03-2020', timePT)):(match('30-06-2021', timePT))]
    fit_timePT <- timePT[(match('07-03-2020', timePT)):(match('30-06-2021', timePT))]
    sample_daysPT <- length(fit_timePT)
    sample_time=1:sample_daysPT
    new_casesPT <- rep(0,sample_daysPT)
    new_casesPT[1] <- cum_confPT[1]
    for (t in 2:sample_daysPT){ 
      new_casesPT[t] <- cum_confPT[t]-cum_confPT[t-1]
    }
    deadPT <- rep(0,sample_daysPT)
    deadPT[1] <- cum_deadPT[1]
    for (t in 2:sample_daysPT){
      deadPT[t] <- cum_deadPT[t]-cum_deadPT[t-1]
      if (deadPT[t]<0) {
        deadPT[t]=0
      }
    }
    vaccines_portugal <- read_csv("https://raw.githubusercontent.com/dssg-pt/covid19pt-data/master/vacinas.csv")
    new_vaccines1 <- vaccines_portugal$doses1_novas
    new_vaccines1[is.na(new_vaccines1)] <- 0
    vaccines_t_45days <- c(rep(0,(length(seq(as.Date('2020/03/07'),as.Date('2020/12/26'),"days")) + 45)),new_vaccines1)
    vaccines_t_45days_PT <- vaccines_t_45days[1:sample_daysPT]
    
    infection_death = rep(0,sample_daysPT)
    infection_death[1] = pgamma(1.5,shape=6.29,rate =0.26) - pgamma(0,shape=6.29,rate =0.26)
    for(i in 2:length(infection_death)) {
      infection_death[i] = pgamma(i+.5,shape=6.29,rate =0.26) - pgamma(i-.5,shape=6.29,rate =0.26)
    }
    initial_PT <- c(10196557,45,45,30,30,0,40)
    left <- seq(1,sample_daysPT,1)
    right <- seq(2,(sample_daysPT+1),1)
    approx_ifr_PT <- c(0.01157665,0.008165473,0.01157665,0.00001,0.000001)
    
    cov_data_by_country = list(n_obs = sample_daysPT,
                    n_difeq = 7,
                    n_pop = population_PT,
                    v = vaccines_t_45days_PT,
                    rho=0.9,
                    yD = deadPT,
                    y_init = initial_PT,
                    t0 = 0,
                    ts = sample_time,
                    I_D = infection_death,
                    left_t=left,
                    right_t=right,
                    time_change_gamma=match('01-01-2021', fit_timePT),
                    sigmaBM_cp1=match('01-08-2020', fit_timePT),
                    sigmaBM_cp2=match('21-12-2020', fit_timePT),
                    ifr_cp1=match('01-06-2020', fit_timePT),
                    ifr_cp2=match('16-11-2020', fit_timePT),  
                    ifr_cp3=match('12-01-2021', fit_timePT),  
                    ifr_cp4=match('01-05-2021', fit_timePT),   
                    ifr_mu=approx_ifr_PT) 
  }
  
  else if (country == "UK"){
    population_UK <- 67886011    #Population numbers: https://www.worldometers.info/world-population/population-by-country/
    time_series_covid19_confirmed_global <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
    time_series_covid19_deaths_global <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv")
    full_cum_confUK <- subset(time_series_covid19_confirmed_global, `Country/Region`=="United Kingdom")[c(5:length(time_series_covid19_confirmed_global))]
    full_cum_deadUK <- subset(time_series_covid19_deaths_global, `Country/Region`=="United Kingdom")[c(5:length(time_series_covid19_deaths_global))]
    aggreg_cum_confUK <- as.matrix(full_cum_confUK)
    aggreg_cum_deadUK <- as.matrix(full_cum_deadUK)
    timeUK <- as.vector(names(full_cum_confUK))
    cum_confUK <- as.numeric(colSums(aggreg_cum_confUK))[(match('2/28/20', timeUK)):(match('6/30/21', timeUK))]
    cum_deadUK <- as.numeric(colSums(aggreg_cum_deadUK))[(match('2/28/20', timeUK)):(match('6/30/21', timeUK))]
    fit_timeUK <- timeUK[(match('2/28/20', timeUK)):(match('6/30/21', timeUK))]
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
      }
    }
    
    vaccines <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv")
    vaccines_UK <- subset(vaccines, `location`=="United Kingdom")
    cum_vaccines1 <- vaccines_UK$people_vaccinated
    cum_vaccines1[is.na(cum_vaccines1)] <- 0
    new_vaccines1 <- rep(0,length(cum_vaccines1))
    new_vaccines1[1:29] <- cum_vaccines1[1:29]
    for (t in 30:length(cum_vaccines1)){
      new_vaccines1[t] <- cum_vaccines1[t]-cum_vaccines1[t-1]
    }
    vaccines_t_45days <- c(rep(0,(length(seq(as.Date('2020/02/28'),as.Date('2020/12/12'),"days")) + 45)),new_vaccines1)
    vaccines_t_45days_UK  <- vaccines_t_45days[1:sample_daysUK]
    
    infection_death <-  rep(0,sample_daysUK)
    infection_death[1] <-  pgamma(1.5,shape=6.29,rate =0.26) - pgamma(0,shape=6.29,rate =0.26)
    for(i in 2:length(infection_death)) {
      infection_death[i] <-  pgamma(i+.5,shape=6.29,rate =0.26) - pgamma(i-.5,shape=6.29,rate =0.26)
    }
    left <- seq(1,sample_daysUK,1)
    right <- seq(2,(sample_daysUK+1),1)
    initial_UK <- c(67884511,450,450,300,300,0,300)
    approx_ifr_UK  <- c(0.01035,0.007245,0.0095,0.00004,0.000018)
    
    cov_data_by_country <- list(n_obs = sample_daysUK,
                                n_difeq = 7, 
                                n_pop = population_UK, 
                                v = vaccines_t_45days_UK , 
                                rho=0.9,
                                yD = deadUK,
                                y_init = initial_UK, 
                                t0 = 0, 
                                ts = sample_time,
                                I_D = infection_death, 
                                left_t=left, 
                                right_t=right, 
                                time_change_gamma=match('1/1/21', fit_timeUK),
                                sigmaBM_cp1=match('8/1/20', fit_timeUK),
                                sigmaBM_cp2=match('11/30/20', fit_timeUK),
                                ifr_cp1=match('7/18/20', fit_timeUK),
                                ifr_cp2=match('10/1/20', fit_timeUK),
                                ifr_cp3=match('1/30/21', fit_timeUK),
                                ifr_cp4=match('6/1/21', fit_timeUK),
                                ifr_mu=approx_ifr_UK) 
  }
  
  else if (country == "DE"){
    population_DE <- 83783942    #Population numbers: https://www.worldometers.info/world-population/population-by-country/
    time_series_covid19_confirmed_global <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
    time_series_covid19_deaths_global <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv")
    full_cum_confDE <- subset(time_series_covid19_confirmed_global, `Country/Region`=="Germany")[c(5:length(time_series_covid19_confirmed_global))]
    full_cum_deadDE <- subset(time_series_covid19_deaths_global, `Country/Region`=="Germany")[c(5:length(time_series_covid19_deaths_global))]
    timeDE <- as.vector(names(full_cum_confDE))
    cum_confDE <- as.numeric(full_cum_confDE)[(match('3/1/20', timeDE)):(match('6/30/21', timeDE))]
    cum_deadDE <- as.numeric(full_cum_deadDE)[(match('3/1/20', timeDE)):(match('6/30/21', timeDE))]
    fit_timeDE <- timeDE[(match('3/1/20', timeDE)):(match('6/30/21', timeDE))]
    sample_daysDE <- length(fit_timeDE)
    sample_time=1:sample_daysDE
    new_casesDE <- rep(0,sample_daysDE)
    new_casesDE[1] <- cum_confDE[1]
    for (t in 2:sample_daysDE){ 
      new_casesDE[t] <- cum_confDE[t]-cum_confDE[t-1]
    }
    deadDE <- rep(0,sample_daysDE)
    deadDE[1] <- cum_deadDE[1]
    for (t in 2:sample_daysDE){
      deadDE[t] <- cum_deadDE[t]-cum_deadDE[t-1]
      if (deadDE[t]<0) {
        deadDE[t]=0
      }
    }
    vaccines <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv")
    vaccines_DE <- subset(vaccines, `location`=="Germany")
    cum_vaccines1 <- vaccines_DE$people_vaccinated
    new_vaccines1 <- rep(0,length(cum_vaccines1))
    new_vaccines1[1] <- cum_vaccines1[1]
    for (t in 2:length(cum_vaccines1)){
      new_vaccines1[t] <- cum_vaccines1[t]-cum_vaccines1[t-1]
    }
    vaccines_t_45days <- c(rep(0,(length(seq(as.Date('2020/03/01'),as.Date('2020/12/26'),"days")) + 45)),new_vaccines1)
    vaccines_t_45days_DE <- vaccines_t_45days[1:sample_daysDE]
    
    infection_death = rep(0,sample_daysDE)
    infection_death[1] = pgamma(1.5,shape=6.29,rate =0.26) - pgamma(0,shape=6.29,rate =0.26)
    for(i in 2:length(infection_death)) {
      infection_death[i] = pgamma(i+.5,shape=6.29,rate =0.26) - pgamma(i-.5,shape=6.29,rate =0.26)
    }
    
    left <- seq(1,sample_daysDE,1)
    right <- seq(2,(sample_daysDE+1),1)
    initial_DE <- c(83782942,300,300,200,200,0,250)
    approx_ifr_DE <- c(0.011039996,0.005434144,0.011420979,0.005,0.0015)
    
    cov_data_by_country = list(n_obs = sample_daysDE,
                    n_difeq = 7,
                    n_pop = population_DE,
                    v = vaccines_t_45days_DE,
                    rho=0.9,
                    yD = deadDE,
                    y_init = initial_DE,
                    t0 = 0,
                    ts = sample_time,
                    I_D = infection_death,
                    left_t=left,
                    right_t=right,
                    time_change_gamma=match('1/1/21', fit_timeDE),
                    sigmaBM_cp1=match('8/1/20', fit_timeDE),
                    sigmaBM_cp2=match('2/20/21', fit_timeDE),
                    ifr_cp1=match('5/27/20', fit_timeDE),
                    ifr_cp2=match('11/25/20', fit_timeDE),
                    ifr_cp3=match('2/24/21', fit_timeDE),
                    ifr_cp4=match('4/16/21', fit_timeDE),
                    ifr_mu=approx_ifr_DE) 
    
  }
  
  else if (country == "SE"){
    population_SE <- 10099265
    # https://www.arcgis.com/sharing/rest/content/items/b5e7488e117749c19881cce45db13f7e/data 
    Sweden_data_C <- read.csv('./data/Sweden_cases.csv', head = TRUE, sep=";")
    Sweden_data_D <- read.csv('./data/Sweden_deaths.csv', head = TRUE, sep=";")
    full_new_casesSE <- Sweden_data_C$Total_Number_of_Cases
    full_deadSE <- Sweden_data_D$Number_of_deaths
    full_deadSE <- c(rep(0,length(seq(as.Date('2020/02/04'),as.Date('2020/03/10'),"days"))),full_deadSE)
    timeSE <- Sweden_data_C$ï.¿date
    fit_timeSE <- timeSE[(match('4/3/2020', timeSE)):(match('30/6/2021', timeSE))]
    new_casesSE <- full_new_casesSE[(match('4/3/2020', timeSE)):(match('30/6/2021', timeSE))]
    deadSE <- full_deadSE[(match('4/3/2020', timeSE)):(match('30/6/2021', timeSE))]
    sample_daysSE <- length(fit_timeSE)
    sample_time=1:sample_daysSE
    
    Sweden_data_V <- read.csv('./data/Sweden_vaccinations.csv', head = TRUE, sep=";") # weekly to daily data
    new_vaccines1 <- Sweden_data_V$daily_vaccinations
    vaccines_t_45days <- c(rep(0,(length(seq(as.Date('2020/03/04'),as.Date('2020/12/26'),"days")) + 45)),new_vaccines1)
    vaccines_t_45days_SE <- vaccines_t_45days[1:sample_daysSE]
  
    infection_death = rep(0,sample_daysSE)
    infection_death[1] = pgamma(1.5,shape=6.29,rate =0.26) - pgamma(0,shape=6.29,rate =0.26)
    for(i in 2:length(infection_death)) {
      infection_death[i] = pgamma(i+.5,shape=6.29,rate =0.26) - pgamma(i-.5,shape=6.29,rate =0.26)
    }
    
    left <- seq(1,sample_daysSE,1)
    right <- seq(2,(sample_daysSE+1),1)
    initial_SE <- c(10098905,100,100,80,80,0,90)
    approx_ifr_SE <- c(0.0103,0.007,0.009,0.001,0.00001)
    
    cov_data_by_country = list(n_obs = sample_daysSE,
                    n_difeq = 7,
                    n_pop = population_SE,
                    v = vaccines_t_45days_SE,
                    rho=0.9,
                    yD = deadSE,
                    y_init = initial_SE,
                    t0 = 0,
                    ts = sample_time,
                    I_D = infection_death,
                    left_t=left,
                    right_t=right,
                    time_change_gamma=match("1/1/2021", fit_timeSE),
                    sigmaBM_cp1=match("1/9/2020", fit_timeSE),
                    sigmaBM_cp2=match("8/2/2021", fit_timeSE),
                    ifr_cp1=match("16/7/2020", fit_timeSE),
                    ifr_cp2=match("2/10/2020", fit_timeSE),
                    ifr_cp3=match("20/12/2020", fit_timeSE),
                    ifr_cp4=match("1/2/2021", fit_timeSE),
                    ifr_mu=approx_ifr_SE) 
    
  }
  
  else if (country == "NO"){
    population_NO <- 5421241
    time_series_covid19_confirmed_global <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
    time_series_covid19_deaths_global <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv")
    full_cum_confNO <- subset(time_series_covid19_confirmed_global, `Country/Region`=="Norway")[c(5:length(time_series_covid19_confirmed_global))]
    full_cum_deadNO <- subset(time_series_covid19_deaths_global, `Country/Region`=="Norway")[c(5:length(time_series_covid19_deaths_global))]
    timeNO <- as.vector(names(full_cum_confNO))
    cum_confNO <- as.numeric(full_cum_confNO)[(match('3/11/20', timeNO)):(match('6/30/21', timeNO))]
    cum_deadNO <- as.numeric(full_cum_deadNO)[(match('3/11/20', timeNO)):(match('6/30/21', timeNO))]
    fit_timeNO <- timeNO[(match('3/11/20', timeNO)):(match('6/30/21', timeNO))]
    sample_daysNO <- length(fit_timeNO)
    sample_time=1:sample_daysNO
    new_casesNO <- rep(0,sample_daysNO)
    new_casesNO[1] <- cum_confNO[1]
    for (t in 2:sample_daysNO){ 
      new_casesNO[t] <- cum_confNO[t]-cum_confNO[t-1]
    }
    deadNO <- rep(0,sample_daysNO)
    deadNO[1] <- cum_deadNO[1]
    for (t in 2:sample_daysNO){
      deadNO[t] <- cum_deadNO[t]-cum_deadNO[t-1]
      if (deadNO[t]<0) {
        deadNO[t]=0
      }
    }
    vaccines <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv")
    vaccines_NO <- subset(vaccines, `location`=="Norway")
    cum_vaccines1 <- vaccines_NO$people_vaccinated
    new_vaccines1 <- rep(0,length(cum_vaccines1))
    new_vaccines1[1] <- cum_vaccines1[1]
    for (t in 2:length(cum_vaccines1)){
      new_vaccines1[t] <- cum_vaccines1[t]-cum_vaccines1[t-1]
    }
    vaccines_t_45days <- c(rep(0,(length(seq(as.Date('2020/03/11'),as.Date('2020/12/01'),"days")) + 45)),new_vaccines1)
    vaccines_t_45days_NO <- vaccines_t_45days[1:sample_daysNO]
    
    infection_death = rep(0,sample_daysNO)
    infection_death[1] = pgamma(1.5,shape=6.29,rate =0.26) - pgamma(0,shape=6.29,rate =0.26)
    for(i in 2:length(infection_death)) {
      infection_death[i] = pgamma(i+.5,shape=6.29,rate =0.26) - pgamma(i-.5,shape=6.29,rate =0.26)
    }
    initial_NO <- c(5420341,250,250,200,200,0,200)
    left <- seq(1,sample_daysNO,1)
    right <- seq(2,(sample_daysNO+1),1)
    approx_ifr_NO <- c(0.0091,0.006,0.005,0.004,0.002)
    
    cov_data_by_country = list(n_obs = sample_daysNO,
                    n_difeq = 7,
                    n_pop = population_NO,
                    v = vaccines_t_45days_NO,
                    rho=0.9,
                    yD = deadNO,
                    y_init = initial_NO,
                    t0 = 0,
                    ts = sample_time,
                    I_D = infection_death,
                    left_t=left,
                    right_t=right,
                    time_change_gamma=match('1/1/21', fit_timeNO),    
                    sigmaBM_cp1=match('7/15/20', fit_timeNO),
                    sigmaBM_cp2=match('2/7/21', fit_timeNO),
                    ifr_cp1=match('5/13/20', fit_timeNO),
                    ifr_cp2=match('10/1/20', fit_timeNO),
                    ifr_cp3=match('12/23/20', fit_timeNO),
                    ifr_cp4=match('1/21/21', fit_timeNO),
                    ifr_mu=approx_ifr_NO) 
    
  }
  
  else{
    print("country not included")
  }

 return(cov_data_by_country) 
}
