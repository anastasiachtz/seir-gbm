reg_data <- function(country = c("GR", "PT", "UK", "DE", "SE", "NO")){
  
  if (country == "GR"){
    load("GR_SEIR_VAC.RData")
    posts_1 <-  rstan::extract(fit_GR)
    post_beta <- posts_1$beta_N
    median_beta = apply(post_beta, 2, median)

    fit_cases <- posts_1$c_tot
    median_fit_cases = apply(fit_cases, 2, median)
    
    time_series_covid19_confirmed_global <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
    full_cum_confGR <- subset(time_series_covid19_confirmed_global, `Country/Region`=="Greece")[c(5:length(time_series_covid19_confirmed_global))]
    timeGR <- as.vector(names(full_cum_confGR))
    cum_confGR <- as.numeric(full_cum_confGR)[(match('3/7/20', timeGR)):(match('9/30/21', timeGR))]
    fit_timeGR <- timeGR[(match('3/7/20', timeGR)):(match('9/30/21', timeGR))]   # 07/03 - 30/09
    sample_daysGR <- length(fit_timeGR)
    new_casesGR <- rep(0,sample_daysGR)
    new_casesGR[1] <- cum_confGR[1]
    for (t in 2:sample_daysGR){ 
      new_casesGR[t] <- cum_confGR[t]-cum_confGR[t-1]
    }
    rescaled_ratio  <- matrix(nrow=2000, ncol=(sample_daysGR))
    for (i in 1:2000){
      for (t in 1:6){
        rescaled_ratio[i,t] <- 0
      }
      for (t in 7:(sample_daysGR)){
        if (fit_cases[i,(t-6)]<new_casesGR[t]) {
          rescaled_ratio[i,t] <- 1
        }
        else
          rescaled_ratio[i,t] <- new_casesGR[t]/fit_cases[i,(t-6)]
      }
    }
    rescaled_ratio <- rescaled_ratio[,7:sample_daysGR]
    median_rescaled_ratio = apply(rescaled_ratio, 2, median)
    time_ratio <- seq(as.Date('2020/03/13'),as.Date('2021/09/30'),"days")
    ppp<-  ggplot(data = data.frame(median_rescaled_ratio, time_ratio), aes(x=time_ratio, y=median_rescaled_ratio)) +
      stat_smooth(method="gam", formula = y ~ s(x), n = length(median_rescaled_ratio))+
      geom_hline(yintercept=1, color = "red")+
      labs(x = "Time", y = "Reported Cases/True Cases") +
      scale_x_date(limits = as.Date(c("2020-03-13","2021-09-30")),labels = date_format("%d-%m"), date_breaks = "14 days") +
      scale_y_continuous(limits=c(0,1), breaks=c(seq(0,1,0.1))) +
      theme_bw(base_size=23)+
      theme(axis.text.x=element_text(angle=90, size=25))+
      theme(axis.text.y=element_text(size=25))+
      theme(legend.position=c(.7, 0.7))+
      theme(legend.title=element_blank())+
      theme(legend.position = "none")
    smooth_fit <- ggplot_build(ppp)$data[[1]]
    smooth_ratio_GAM <- smooth_fit$y
    
    ## 19/03/20 - 31/03/21
    log_b <- log(median_beta[(match('3/19/20', fit_timeGR)):(match('3/31/21', fit_timeGR))])
    estimated_total_cases <- median_fit_cases[(match('3/19/20', fit_timeGR)):(match('3/31/21', fit_timeGR))]
    logit_smooth_ratio <- logit(smooth_ratio_GAM[(match(as.Date("2020-03-19"), time_ratio)):(match(as.Date("2021-03-31"), time_ratio))])
    
    dependent <- cbind(estimated_total_cases, log_b, logit_smooth_ratio)
    
    #############################################################################################################
    # Covariates
    
    mobility <- read.csv('./data/Global_Mobility_Report.csv', head = TRUE, sep=",")  # https://www.google.com/covid19/mobility/
    GR_mobility <- subset(mobility,country_region=="Greece")
    aggregate_mobility <- subset(GR_mobility, sub_region_1=="")
    time_mobility <- aggregate_mobility$date
    mob_pca <- prcomp(aggregate_mobility[c(10:15)], center = TRUE, scale = TRUE)
    summary(mob_pca)
    firstPCmobAll <-  (mob_pca$x[,1])*(-1)
    
    # Discretize serial interval distribution
    serial_int <- rep(0,length(firstPCmobAll))
    serial_int_rev <- rep(0,length(firstPCmobAll))
    weighted_firstPCmobAll <- rep(0,length(firstPCmobAll))

    serial_int[1] <- pgamma(1.5,shape=2.6,rate =0.4) - pgamma(0,shape=2.6,rate =0.4)
    for(i in 2:length(serial_int)) {
      serial_int[i] <- pgamma(i+.5,shape=2.6,rate =0.4) - pgamma(i-.5,shape=2.6,rate =0.4)
    }
    for (i in 1:length(firstPCmobAll)){
        serial_int_rev[i] <- serial_int[length(firstPCmobAll)-i+1]
    }
    for (i in 2:length(firstPCmobAll)){
      weighted_firstPCmobAll[i] <- (head(firstPCmobAll,i-1)) %*% (tail(serial_int_rev,i-1))
    }
    weighted_mob <- weighted_firstPCmobAll[(match("2020-03-19", time_mobility)):(match("2021-03-31", time_mobility))]
    
    covid_tests <- read.csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/testing/covid-testing-all-observations.csv")
    GR <- subset(covid_tests, `Entity`=="Greece - samples tested")
    test_timeGR <- as.Date(GR$Date)
    tests_GR <- GR$Daily.change.in.cumulative.total
    for (i in 11:length(tests_GR)){
      if (is.na(tests_GR[i])){
        tests_GR[i] <- tests_GR[i-1]
      }
    }
    test_t_3 <- tests_GR[(match(as.Date("2020-03-16"), test_timeGR)):(match(as.Date("2021-03-28"), test_timeGR))]
    test_t_4 <- tests_GR[(match(as.Date("2020-03-15"), test_timeGR)):(match(as.Date("2021-03-27"), test_timeGR))]
    test_t_5 <- tests_GR[(match(as.Date("2020-03-14"), test_timeGR)):(match(as.Date("2021-03-26"), test_timeGR))]
    test_t_6 <- tests_GR[(match(as.Date("2020-03-13"), test_timeGR)):(match(as.Date("2021-03-25"), test_timeGR))]
    tests <- cbind(test_t_3, test_t_4, test_t_5, test_t_6)
    
    covariates <- cbind(weighted_mob, tests)
    XX <- t(covariates)%*%covariates
    inv_XX <- chol2inv(chol(XX))
    mean_g <- rep(0, ncol(covariates))
    
    # Modify data into a form suitable for Stan
    reg_data_by_country = list(n_obs = nrow(dependent),
                    K = ncol(dependent),  #number of regressions
                    J = ncol(covariates), #number of covariates
                    x = covariates,
                    y = dependent,
                    b0 = mean_g,
                    inv_XX = inv_XX)
  }
  
  else if (country == "PT"){
    load("PT_SEIR_VAC.RData")
    posts_1 <-  rstan::extract(fit_PT)
    post_beta <- posts_1$beta_N
    median_beta = apply(post_beta, 2, median)
    
    fit_cases <- posts_1$c_tot
    median_fit_cases = apply(fit_cases, 2, median)
    
    portugal_data <- read_csv("https://raw.githubusercontent.com/dssg-pt/covid19pt-data/master/data.csv")
    full_cum_confPT <- portugal_data$confirmados
    timePT <- portugal_data$data
    cum_confPT <- as.numeric(full_cum_confPT)[(match('07-03-2020', timePT)):(match('30-09-2021', timePT))]
    fit_timePT <- timePT[(match('07-03-2020', timePT)):(match('30-09-2021', timePT))]
    sample_daysPT <- length(fit_timePT)
    new_casesPT <- rep(0,sample_daysPT)
    new_casesPT[1] <- cum_confPT[1]
    for (t in 2:sample_daysPT){ 
      new_casesPT[t] <- cum_confPT[t]-cum_confPT[t-1]
    }
    rescaled_ratio  <- matrix(nrow=2000, ncol=(sample_daysPT))
    for (i in 1:2000){
      for (t in 1:6){
        rescaled_ratio[i,t] <- 0
      }
      for (t in 7:(sample_daysPT)){
        if (fit_cases[i,(t-6)]<new_casesPT[t]) {
          rescaled_ratio[i,t] <- 1
        }
        else
          rescaled_ratio[i,t] <- new_casesPT[t]/fit_cases[i,(t-6)]
      }
    }
    rescaled_ratio <- rescaled_ratio[,7:sample_daysPT]
    median_rescaled_ratio = apply(rescaled_ratio, 2, median)
    time_ratio <- seq(as.Date('2020/03/13'),as.Date('2021/09/30'),"days")
    ppp<-  ggplot(data = data.frame(median_rescaled_ratio, time_ratio), aes(x=time_ratio, y=median_rescaled_ratio)) +
      stat_smooth(method="gam", formula = y ~ s(x), n = length(median_rescaled_ratio))+
      geom_hline(yintercept=1, color = "red")+
      labs(x = "Time", y = "Reported Cases/True Cases") +
      scale_x_date(limits = as.Date(c("2020-03-13","2021-09-30")),labels = date_format("%d-%m"), date_breaks = "14 days") +
      scale_y_continuous(limits=c(0,1), breaks=c(seq(0,1,0.1))) +
      theme_bw(base_size=23)+
      theme(axis.text.x=element_text(angle=90, size=25))+
      theme(axis.text.y=element_text(size=25))+
      theme(legend.position=c(.7, 0.7))+
      theme(legend.title=element_blank())+
      theme(legend.position = "none")
    smooth_fit <- ggplot_build(ppp)$data[[1]]
    smooth_ratio_GAM <- smooth_fit$y
    
    ## 13/03/20 - 31/03/21
    log_b <- log(median_beta[(match('13-03-2020', fit_timePT)):(match('31-03-2021', fit_timePT))])
    estimated_total_cases <- median_fit_cases[(match('13-03-2020', fit_timePT)):(match('31-03-2021', fit_timePT))]
    logit_smooth_ratio <- logit(smooth_ratio_GAM[(match(as.Date("2020-03-13"), time_ratio)):(match(as.Date("2021-03-31"), time_ratio))])
    
    dependent <- cbind(estimated_total_cases, log_b, logit_smooth_ratio)
    
    #############################################################################################################
    # Covariates
    
    mobility <- read.csv('./data/Global_Mobility_Report.csv', head = TRUE, sep=",")  # https://www.google.com/covid19/mobility/
    PT_mobility <- subset(mobility,country_region=="Portugal")
    aggregate_mobility <- subset(PT_mobility, sub_region_1=="")
    time_mobility <- aggregate_mobility$date
    mob_pca <- prcomp(aggregate_mobility[c(10:15)], center = TRUE, scale = TRUE)
    summary(mob_pca)
    firstPCmobAll <-  (mob_pca$x[,1])*(-1)
    
    # Discretize serial interval distribution
    serial_int <- rep(0,length(firstPCmobAll))
    serial_int_rev <- rep(0,length(firstPCmobAll))
    weighted_firstPCmobAll <- rep(0,length(firstPCmobAll))
    
    serial_int[1] <- pgamma(1.5,shape=2.6,rate =0.4) - pgamma(0,shape=2.6,rate =0.4)
    for(i in 2:length(serial_int)) {
      serial_int[i] <- pgamma(i+.5,shape=2.6,rate =0.4) - pgamma(i-.5,shape=2.6,rate =0.4)
    }
    for (i in 1:length(firstPCmobAll)){
      serial_int_rev[i] <- serial_int[length(firstPCmobAll)-i+1]
    }
    for (i in 2:length(firstPCmobAll)){
      weighted_firstPCmobAll[i] <- (head(firstPCmobAll,i-1)) %*% (tail(serial_int_rev,i-1))
    }
    weighted_mob <- weighted_firstPCmobAll[(match("2020-03-13", time_mobility)):(match("2021-03-31", time_mobility))]
    
    covid_tests <- read.csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/testing/covid-testing-all-observations.csv")
    PT <- subset(covid_tests, `Entity`=="Portugal - tests performed")
    test_timePT <- as.Date(PT$Date)
    tests_PT <- PT$Daily.change.in.cumulative.total
    test_t_3 <- tests_PT[(match(as.Date("2020-03-10"), test_timePT)):(match(as.Date("2021-03-28"), test_timePT))]
    test_t_4 <- tests_PT[(match(as.Date("2020-03-09"), test_timePT)):(match(as.Date("2021-03-27"), test_timePT))]
    test_t_5 <- tests_PT[(match(as.Date("2020-03-08"), test_timePT)):(match(as.Date("2021-03-26"), test_timePT))]
    test_t_6 <- tests_PT[(match(as.Date("2020-03-07"), test_timePT)):(match(as.Date("2021-03-25"), test_timePT))]
    tests <- cbind(test_t_3, test_t_4, test_t_5, test_t_6)
    
    covariates <- cbind(weighted_mob, tests)
    XX <- t(covariates)%*%covariates
    inv_XX <- chol2inv(chol(XX))
    mean_g <- rep(0, ncol(covariates))
    
    # Modify data into a form suitable for Stan
    reg_data_by_country = list(n_obs = nrow(dependent),
                    K = ncol(dependent),  #number of regressions
                    J = ncol(covariates), #number of covariates
                    x = covariates,
                    y = dependent,
                    b0 = mean_g,
                    inv_XX = inv_XX) 
  }
  
  else if (country == "UK"){
    load("UK_SEIR_VAC.RData")
    posts_1 <-  rstan::extract(fit_UK)
    post_beta <- posts_1$beta_N
    median_beta = apply(post_beta, 2, median)
    
    fit_cases <- posts_1$c_tot
    median_fit_cases = apply(fit_cases, 2, median)
    
    time_series_covid19_confirmed_global <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
    full_cum_confUK <- subset(time_series_covid19_confirmed_global, `Country/Region`=="United Kingdom")[c(5:length(time_series_covid19_confirmed_global))]
    aggreg_cum_confUK <- as.matrix(full_cum_confUK)
    timeUK <- as.vector(names(full_cum_confUK))
    cum_confUK <- as.numeric(colSums(aggreg_cum_confUK))[(match('2/28/20', timeUK)):(match('9/30/21', timeUK))]
    fit_timeUK <- timeUK[(match('2/28/20', timeUK)):(match('9/30/21', timeUK))]
    sample_daysUK <- length(fit_timeUK)
    new_casesUK <- rep(0,sample_daysUK)
    new_casesUK[1] <- cum_confUK[1]
    for (t in 2:sample_daysUK){ 
      new_casesUK[t] <- cum_confUK[t]-cum_confUK[t-1]
      if (new_casesUK[t]<0) {
        new_casesUK[t]<-new_casesUK[t-1]
      }
    }
    rescaled_ratio  <- matrix(nrow=2000, ncol=(sample_daysUK))
    for (i in 1:2000){
      for (t in 1:6){
        rescaled_ratio[i,t] <- 0
      }
      for (t in 7:(sample_daysUK)){
        if (fit_cases[i,(t-6)]<new_casesUK[t]) {
          rescaled_ratio[i,t] <- 1
        }
        else
          rescaled_ratio[i,t] <- new_casesUK[t]/fit_cases[i,(t-6)]
      }
    }
    rescaled_ratio <- rescaled_ratio[,7:sample_daysUK]
    median_rescaled_ratio = apply(rescaled_ratio, 2, median)
    time_ratio <- seq(as.Date('2020/03/05'),as.Date('2021/09/30'),"days")
    ppp<-  ggplot(data = data.frame(median_rescaled_ratio, time_ratio), aes(x=time_ratio, y=median_rescaled_ratio)) +
      stat_smooth(method="gam", formula = y ~ s(x), n = length(median_rescaled_ratio))+
      geom_hline(yintercept=1, color = "red")+
      labs(x = "Time", y = "Reported Cases/True Cases") +
      scale_x_date(limits = as.Date(c("2020-03-05","2021-09-30")),labels = date_format("%d-%m"), date_breaks = "14 days") +
      scale_y_continuous(limits=c(0,1), breaks=c(seq(0,1,0.1))) +
      theme_bw(base_size=23)+
      theme(axis.text.x=element_text(angle=90, size=25))+
      theme(axis.text.y=element_text(size=25))+
      theme(legend.position=c(.7, 0.7))+
      theme(legend.title=element_blank())+
      theme(legend.position = "none")
    smooth_fit <- ggplot_build(ppp)$data[[1]]
    smooth_ratio_GAM <- smooth_fit$y
    
    ## 07/04/20 - 31/03/21
    log_b <- log(median_beta[(match('4/7/20', fit_timeUK)):(match('3/31/21', fit_timeUK))])
    estimated_total_cases <- median_fit_cases[(match('4/7/20', fit_timeUK)):(match('3/31/21', fit_timeUK))]
    logit_smooth_ratio <- logit(smooth_ratio_GAM[(match(as.Date("2020-04-07"), time_ratio)):(match(as.Date("2021-03-31"), time_ratio))])
    
    dependent <- cbind(estimated_total_cases, log_b, logit_smooth_ratio)
    
    #############################################################################################################
    # Covariates
    
    mobility <- read.csv('./data/Global_Mobility_Report.csv', head = TRUE, sep=",")  # https://www.google.com/covid19/mobility/
    UK_mobility <- subset(mobility,country_region=="United Kingdom")
    aggregate_mobility <- subset(UK_mobility, sub_region_1=="")
    time_mobility <- aggregate_mobility$date
    mob_pca <- prcomp(aggregate_mobility[c(10:15)], center = TRUE, scale = TRUE)
    summary(mob_pca)
    firstPCmobAll <-  (mob_pca$x[,1])*(-1)
    
    # Discretize serial interval distribution
    serial_int <- rep(0,length(firstPCmobAll))
    serial_int_rev <- rep(0,length(firstPCmobAll))
    weighted_firstPCmobAll <- rep(0,length(firstPCmobAll))
    
    serial_int[1] <- pgamma(1.5,shape=2.6,rate =0.4) - pgamma(0,shape=2.6,rate =0.4)
    for(i in 2:length(serial_int)) {
      serial_int[i] <- pgamma(i+.5,shape=2.6,rate =0.4) - pgamma(i-.5,shape=2.6,rate =0.4)
    }
    for (i in 1:length(firstPCmobAll)){
      serial_int_rev[i] <- serial_int[length(firstPCmobAll)-i+1]
    }
    for (i in 2:length(firstPCmobAll)){
      weighted_firstPCmobAll[i] <- (head(firstPCmobAll,i-1)) %*% (tail(serial_int_rev,i-1))
    }
    weighted_mob <- weighted_firstPCmobAll[(match("2020-04-07", time_mobility)):(match("2021-03-31", time_mobility))]
    
    covid_tests <- read.csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/testing/covid-testing-all-observations.csv")
    UK <- subset(covid_tests, `Entity`=="United Kingdom - tests performed")
    test_timeUK <- as.Date(UK$Date)
    tests_UK <- UK$Daily.change.in.cumulative.total
    test_t_3 <- tests_UK[(match(as.Date("2020-04-04"), test_timeUK)):(match(as.Date("2021-03-28"), test_timeUK))]
    test_t_4 <- tests_UK[(match(as.Date("2020-04-03"), test_timeUK)):(match(as.Date("2021-03-27"), test_timeUK))]
    test_t_5 <- tests_UK[(match(as.Date("2020-04-02"), test_timeUK)):(match(as.Date("2021-03-26"), test_timeUK))]
    test_t_6 <- tests_UK[(match(as.Date("2020-04-01"), test_timeUK)):(match(as.Date("2021-03-25"), test_timeUK))]
    tests <- cbind(test_t_3, test_t_4, test_t_5, test_t_6)
    
    covariates <- cbind(weighted_mob, tests)
    XX <- t(covariates)%*%covariates
    inv_XX <- chol2inv(chol(XX))
    mean_g <- rep(0, ncol(covariates))
    
    # Modify data into a form suitable for Stan
    reg_data_by_country = list(n_obs = nrow(dependent),
                               K = ncol(dependent),  #number of regressions
                               J = ncol(covariates), #number of covariates
                               x = covariates,
                               y = dependent,
                               b0 = mean_g,
                               inv_XX = inv_XX)
  }
  
  else if (country == "DE"){
    load("DE_SEIR_VAC.RData")
    posts_1 <-  rstan::extract(fit_DE)
    post_beta <- posts_1$beta_N
    median_beta = apply(post_beta, 2, median)
    
    fit_cases <- posts_1$c_tot
    median_fit_cases = apply(fit_cases, 2, median)
    
    time_series_covid19_confirmed_global <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
    full_cum_confDE <- subset(time_series_covid19_confirmed_global, `Country/Region`=="Germany")[c(5:length(time_series_covid19_confirmed_global))]
    timeDE <- as.vector(names(full_cum_confDE))
    cum_confDE <- as.numeric(full_cum_confDE)[(match('3/1/20', timeDE)):(match('9/30/21', timeDE))]
    fit_timeDE <- timeDE[(match('3/1/20', timeDE)):(match('9/30/21', timeDE))]
    sample_daysDE <- length(fit_timeDE)
    new_casesDE <- rep(0,sample_daysDE)
    new_casesDE[1] <- cum_confDE[1]
    for (t in 2:sample_daysDE){ 
      new_casesDE[t] <- cum_confDE[t]-cum_confDE[t-1]
    }
    rescaled_ratio  <- matrix(nrow=2000, ncol=(sample_daysDE))
    for (i in 1:2000){
      for (t in 1:6){
        rescaled_ratio[i,t] <- 0
      }
      for (t in 7:(sample_daysDE)){
        if (fit_cases[i,(t-6)]<new_casesDE[t]) {
          rescaled_ratio[i,t] <- 1
        }
        else
          rescaled_ratio[i,t] <- new_casesDE[t]/fit_cases[i,(t-6)]
      }
    }
    rescaled_ratio <- rescaled_ratio[,7:sample_daysDE]
    median_rescaled_ratio = apply(rescaled_ratio, 2, median)
    time_ratio <- seq(as.Date('2020/03/07'),as.Date('2021/09/30'),"days")
    ppp<-  ggplot(data = data.frame(median_rescaled_ratio, time_ratio), aes(x=time_ratio, y=median_rescaled_ratio)) +
      stat_smooth(method="gam", formula = y ~ s(x), n = length(median_rescaled_ratio))+
      geom_hline(yintercept=1, color = "red")+
      labs(x = "Time", y = "Reported Cases/True Cases") +
      scale_x_date(limits = as.Date(c("2020-03-07","2021-09-30")),labels = date_format("%d-%m"), date_breaks = "14 days") +
      scale_y_continuous(limits=c(0,1), breaks=c(seq(0,1,0.1))) +
      theme_bw(base_size=23)+
      theme(axis.text.x=element_text(angle=90, size=25))+
      theme(axis.text.y=element_text(size=25))+
      theme(legend.position=c(.7, 0.7))+
      theme(legend.title=element_blank())+
      theme(legend.position = "none")
    smooth_fit <- ggplot_build(ppp)$data[[1]]
    smooth_ratio_GAM <- smooth_fit$y
    
    ## 08/03/20 - 31/03/21
    log_b <- log(median_beta[(match('3/8/20', fit_timeDE)):(match('3/31/21', fit_timeDE))])
    estimated_total_cases <- median_fit_cases[(match('3/8/20', fit_timeDE)):(match('3/31/21', fit_timeDE))]
    logit_smooth_ratio <- logit(smooth_ratio_GAM[(match(as.Date("2020-03-08"), time_ratio)):(match(as.Date("2021-03-31"), time_ratio))])
    
    dependent <- cbind(estimated_total_cases, log_b, logit_smooth_ratio)
    
    #############################################################################################################
    # Covariates
    
    mobility <- read.csv('./data/Global_Mobility_Report.csv', head = TRUE, sep=",")  # https://www.google.com/covid19/mobility/
    DE_mobility <- subset(mobility,country_region=="Germany")
    aggregate_mobility <- subset(DE_mobility, sub_region_1=="")
    time_mobility <- aggregate_mobility$date
    mob_pca <- prcomp(aggregate_mobility[c(10:15)], center = TRUE, scale = TRUE)
    summary(mob_pca)
    firstPCmobAll <-  (mob_pca$x[,1])*(-1)
    
    # Discretize serial interval distribution
    serial_int <- rep(0,length(firstPCmobAll))
    serial_int_rev <- rep(0,length(firstPCmobAll))
    weighted_firstPCmobAll <- rep(0,length(firstPCmobAll))
    
    serial_int[1] <- pgamma(1.5,shape=2.6,rate =0.4) - pgamma(0,shape=2.6,rate =0.4)
    for(i in 2:length(serial_int)) {
      serial_int[i] <- pgamma(i+.5,shape=2.6,rate =0.4) - pgamma(i-.5,shape=2.6,rate =0.4)
    }
    for (i in 1:length(firstPCmobAll)){
      serial_int_rev[i] <- serial_int[length(firstPCmobAll)-i+1]
    }
    for (i in 2:length(firstPCmobAll)){
      weighted_firstPCmobAll[i] <- (head(firstPCmobAll,i-1)) %*% (tail(serial_int_rev,i-1))
    }
    weighted_mob <- weighted_firstPCmobAll[(match("2020-03-08", time_mobility)):(match("2021-03-31", time_mobility))]
    
    DE <- read.csv('./data/Germany_tests.csv', head = TRUE, sep=";")  #https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Testzahl.html;jsessionid=821BB52428F222307CBE3BA46A8A4106.internet091?nn=13490888
    tests_weekly_DE <- DE$Number.of.tests
    test_timeDE <- seq(as.Date('2020/03/02'),as.Date('2021/06/27'),'days')
    tests_DE <- rep(0,length(test_timeDE))
    for (i in 1:length(test_timeDE)) {
      tests_DE[i] <- tests_weekly_DE[1+i/8]/7
      }
    test_t_3 <- tests_DE[(match(as.Date("2020-03-05"), test_timeDE)):(match(as.Date("2021-03-28"), test_timeDE))]
    test_t_4 <- tests_DE[(match(as.Date("2020-03-04"), test_timeDE)):(match(as.Date("2021-03-27"), test_timeDE))]
    test_t_5 <- tests_DE[(match(as.Date("2020-03-03"), test_timeDE)):(match(as.Date("2021-03-26"), test_timeDE))]
    test_t_6 <- tests_DE[(match(as.Date("2020-03-02"), test_timeDE)):(match(as.Date("2021-03-25"), test_timeDE))]
    tests <- cbind(test_t_3, test_t_4, test_t_5, test_t_6)
    
    covariates <- cbind(weighted_mob, tests)
    XX <- t(covariates)%*%covariates
    inv_XX <- chol2inv(chol(XX))
    mean_g <- rep(0, ncol(covariates))
    
    # Modify data into a form suitable for Stan
    reg_data_by_country = list(n_obs = nrow(dependent),
                               K = ncol(dependent),  #number of regressions
                               J = ncol(covariates), #number of covariates
                               x = covariates,
                               y = dependent,
                               b0 = mean_g,
                               inv_XX = inv_XX) 
  }
  
  else if (country == "SE"){
    load("SE_SEIR_VAC.RData")
    posts_1 <-  rstan::extract(fit_SE)
    post_beta <- posts_1$beta_N
    median_beta = apply(post_beta, 2, median)
    
    fit_cases <- posts_1$c_tot
    median_fit_cases = apply(fit_cases, 2, median)
    
    Sweden_data_C <- read.csv('./data/Sweden_cases.csv', head = TRUE, sep=";")
    full_new_casesSE <- Sweden_data_C$Total_Number_of_Cases
    timeSE <- Sweden_data_C$date
    fit_timeSE <- timeSE[(match('4/3/20', timeSE)):(match('30/9/21', timeSE))]
    new_casesSE <- full_new_casesSE[(match('4/3/20', timeSE)):(match('30/9/21', timeSE))]
    sample_daysSE <- length(fit_timeSE)
    rescaled_ratio  <- matrix(nrow=2000, ncol=(sample_daysSE))
    for (i in 1:2000){
      for (t in 1:6){
        rescaled_ratio[i,t] <- 0
      }
      for (t in 7:(sample_daysSE)){
        if (fit_cases[i,(t-6)]<new_casesSE[t]) {
          rescaled_ratio[i,t] <- 1
        }
        else
          rescaled_ratio[i,t] <- new_casesSE[t]/fit_cases[i,(t-6)]
      }
    }
    rescaled_ratio <- rescaled_ratio[,7:sample_daysSE]
    median_rescaled_ratio = apply(rescaled_ratio, 2, median)
    time_ratio <- seq(as.Date('2020/03/10'),as.Date('2021/09/30'),"days")
    ppp<-  ggplot(data = data.frame(median_rescaled_ratio, time_ratio), aes(x=time_ratio, y=median_rescaled_ratio)) +
      stat_smooth(method="gam", formula = y ~ s(x), n = length(median_rescaled_ratio))+
      geom_hline(yintercept=1, color = "red")+
      labs(x = "Time", y = "Reported Cases/True Cases") +
      scale_x_date(limits = as.Date(c("2020-03-10","2021-09-30")),labels = date_format("%d-%m"), date_breaks = "14 days") +
      scale_y_continuous(limits=c(0,1), breaks=c(seq(0,1,0.1))) +
      theme_bw(base_size=23)+
      theme(axis.text.x=element_text(angle=90, size=25))+
      theme(axis.text.y=element_text(size=25))+
      theme(legend.position=c(.7, 0.7))+
      theme(legend.title=element_blank())+
      theme(legend.position = "none")
    smooth_fit <- ggplot_build(ppp)$data[[1]]
    smooth_ratio_GAM <- smooth_fit$y
    for (i in 1:length(smooth_ratio_GAM)){
      if (is.na(smooth_ratio_GAM[i])){
        smooth_ratio_GAM[i] <- 0.9999
      }
    }
    
    ## 19/07/20 - 31/03/21
    log_b <- log(median_beta[(match('19/7/20', fit_timeSE)):(match('31/3/21', fit_timeSE))])
    estimated_total_cases <- median_fit_cases[(match('19/7/20', fit_timeSE)):(match('31/3/21', fit_timeSE))]
    logit_smooth_ratio <- logit(smooth_ratio_GAM[(match(as.Date("2020-07-19"), time_ratio)):(match(as.Date("2021-03-31"), time_ratio))])
    
    dependent <- cbind(estimated_total_cases, log_b, logit_smooth_ratio)
    
    #############################################################################################################
    # Covariates
    
    mobility <- read.csv('./data/Global_Mobility_Report.csv', head = TRUE, sep=",")  # https://www.google.com/covid19/mobility/
    SE_mobility <- subset(mobility,country_region=="Sweden")
    aggregate_mobility <- subset(SE_mobility, sub_region_1=="")
    time_mobility <- aggregate_mobility$date
    mob_pca <- prcomp(aggregate_mobility[c(10:15)], center = TRUE, scale = TRUE)
    summary(mob_pca)
    firstPCmobAll <-  (mob_pca$x[,1])*(-1)
    
    # Discretize serial interval distribution
    serial_int <- rep(0,length(firstPCmobAll))
    serial_int_rev <- rep(0,length(firstPCmobAll))
    weighted_firstPCmobAll <- rep(0,length(firstPCmobAll))
    
    serial_int[1] <- pgamma(1.5,shape=2.6,rate =0.4) - pgamma(0,shape=2.6,rate =0.4)
    for(i in 2:length(serial_int)) {
      serial_int[i] <- pgamma(i+.5,shape=2.6,rate =0.4) - pgamma(i-.5,shape=2.6,rate =0.4)
    }
    for (i in 1:length(firstPCmobAll)){
      serial_int_rev[i] <- serial_int[length(firstPCmobAll)-i+1]
    }
    for (i in 2:length(firstPCmobAll)){
      weighted_firstPCmobAll[i] <- (head(firstPCmobAll,i-1)) %*% (tail(serial_int_rev,i-1))
    }
    weighted_mob <- weighted_firstPCmobAll[(match("2020-07-19", time_mobility)):(match("2021-03-31", time_mobility))]
    
    covid_tests <- read.csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/testing/covid-testing-all-observations.csv")
    SE <- subset(covid_tests, `Entity`=="Sweden - tests performed")
    test_timeSE <- as.Date(SE$Date)
    tests_SE <- SE$Daily.change.in.cumulative.total
    test_t_3 <- tests_SE[(match(as.Date("2020-07-16"), test_timeSE)):(match(as.Date("2021-03-28"), test_timeSE))]
    test_t_4 <- tests_SE[(match(as.Date("2020-07-15"), test_timeSE)):(match(as.Date("2021-03-27"), test_timeSE))]
    test_t_5 <- tests_SE[(match(as.Date("2020-07-14"), test_timeSE)):(match(as.Date("2021-03-26"), test_timeSE))]
    test_t_6 <- tests_SE[(match(as.Date("2020-07-13"), test_timeSE)):(match(as.Date("2021-03-25"), test_timeSE))]
    tests <- cbind(test_t_3, test_t_4, test_t_5, test_t_6)
    
    covariates <- cbind(weighted_mob, tests)
    XX <- t(covariates)%*%covariates
    inv_XX <- chol2inv(chol(XX))
    mean_g <- rep(0, ncol(covariates))
    
    # Modify data into a form suitable for Stan
    reg_data_by_country = list(n_obs = nrow(dependent),
                               K = ncol(dependent),  #number of regressions
                               J = ncol(covariates), #number of covariates
                               x = covariates,
                               y = dependent,
                               b0 = mean_g,
                               inv_XX = inv_XX)
  }
  
  else if (country == "NO"){
    load("NO_SEIR_VAC.RData")
    posts_1 <-  rstan::extract(fit_NO)
    post_beta <- posts_1$beta_N
    median_beta = apply(post_beta, 2, median)
    
    fit_cases <- posts_1$c_tot
    median_fit_cases = apply(fit_cases, 2, median)
    
    time_series_covid19_confirmed_global <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
    full_cum_confNO <- subset(time_series_covid19_confirmed_global, `Country/Region`=="Norway")[c(5:length(time_series_covid19_confirmed_global))]
    timeNO <- as.vector(names(full_cum_confNO))
    cum_confNO <- as.numeric(full_cum_confNO)[(match('3/11/20', timeNO)):(match('9/30/21', timeNO))]
    fit_timeNO <- timeNO[(match('3/11/20', timeNO)):(match('9/30/21', timeNO))]
    sample_daysNO <- length(fit_timeNO)
    new_casesNO <- rep(0,sample_daysNO)
    new_casesNO[1] <- cum_confNO[1]
    for (t in 2:sample_daysNO){ 
      new_casesNO[t] <- cum_confNO[t]-cum_confNO[t-1]
    }
    rescaled_ratio  <- matrix(nrow=2000, ncol=(sample_daysNO))
    for (i in 1:2000){
      for (t in 1:6){
        rescaled_ratio[i,t] <- 0
      }
      for (t in 7:(sample_daysNO)){
        if (fit_cases[i,(t-6)]<new_casesNO[t]) {
          rescaled_ratio[i,t] <- 1
        }
        else
          rescaled_ratio[i,t] <- new_casesNO[t]/fit_cases[i,(t-6)]
      }
    }
    rescaled_ratio <- rescaled_ratio[,7:sample_daysNO]
    median_rescaled_ratio = apply(rescaled_ratio, 2, median)
    time_ratio <- seq(as.Date('2020/03/17'),as.Date('2021/09/30'),"days")
    ppp<-  ggplot(data = data.frame(median_rescaled_ratio, time_ratio), aes(x=time_ratio, y=median_rescaled_ratio)) +
      stat_smooth(method="gam", formula = y ~ s(x), n = length(median_rescaled_ratio))+
      geom_hline(yintercept=1, color = "red")+
      labs(x = "Time", y = "Reported Cases/True Cases") +
      scale_x_date(limits = as.Date(c("2020-03-17","2021-09-30")),labels = date_format("%d-%m"), date_breaks = "14 days") +
      scale_y_continuous(limits=c(0,1), breaks=c(seq(0,1,0.1))) +
      theme_bw(base_size=23)+
      theme(axis.text.x=element_text(angle=90, size=25))+
      theme(axis.text.y=element_text(size=25))+
      theme(legend.position=c(.7, 0.7))+
      theme(legend.title=element_blank())+
      theme(legend.position = "none")
    smooth_fit <- ggplot_build(ppp)$data[[1]]
    smooth_ratio_GAM <- smooth_fit$y
    for (i in 1:length(smooth_ratio_GAM)){
      if (is.na(smooth_ratio_GAM[i])){
        smooth_ratio_GAM[i] <- 0.9999
      }
    }
    
    ## 07/04/20 - 31/03/21
    log_b <- log(median_beta[(match('4/7/20', fit_timeNO)):(match('3/31/21', fit_timeNO))])
    estimated_total_cases <- median_fit_cases[(match('4/7/20', fit_timeNO)):(match('3/31/21', fit_timeNO))]
    logit_smooth_ratio <- logit(smooth_ratio_GAM[(match(as.Date("2020-04-07"), time_ratio)):(match(as.Date("2021-03-31"), time_ratio))])
    
    dependent <- cbind(estimated_total_cases, log_b, logit_smooth_ratio)
    
    #############################################################################################################
    # Covariates
    
    mobility <- read.csv('./data/Global_Mobility_Report.csv', head = TRUE, sep=",")  # https://www.google.com/covid19/mobility/
    NO_mobility <- subset(mobility,country_region=="Norway")
    aggregate_mobility <- subset(NO_mobility, sub_region_1=="")
    time_mobility <- aggregate_mobility$date
    mod_aggregate_mobility <- as.data.frame(aggregate_mobility[c(10:15)])
    mod_aggregate_mobility[is.na(mod_aggregate_mobility)] <- 0
    mob_pca <- prcomp(mod_aggregate_mobility, center = TRUE, scale = TRUE)
    summary(mob_pca)
    firstPCmobAll <-  (mob_pca$x[,1])*(-1)
    
    # Discretize serial interval distribution
    serial_int <- rep(0,length(firstPCmobAll))
    serial_int_rev <- rep(0,length(firstPCmobAll))
    weighted_firstPCmobAll <- rep(0,length(firstPCmobAll))
    
    serial_int[1] <- pgamma(1.5,shape=2.6,rate =0.4) - pgamma(0,shape=2.6,rate =0.4)
    for(i in 2:length(serial_int)) {
      serial_int[i] <- pgamma(i+.5,shape=2.6,rate =0.4) - pgamma(i-.5,shape=2.6,rate =0.4)
    }
    for (i in 1:length(firstPCmobAll)){
      serial_int_rev[i] <- serial_int[length(firstPCmobAll)-i+1]
    }
    for (i in 2:length(firstPCmobAll)){
      weighted_firstPCmobAll[i] <- (head(firstPCmobAll,i-1)) %*% (tail(serial_int_rev,i-1))
    }
    weighted_mob <- weighted_firstPCmobAll[(match("2020-04-07", time_mobility)):(match("2021-03-31", time_mobility))]
    
    covid_tests <- read.csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/testing/covid-testing-all-observations.csv")
    NO <- subset(covid_tests, `Entity`=="Norway - people tested")
    test_timeNO <- as.Date(NO$Date)
    tests_NO <- NO$Daily.change.in.cumulative.total
    for (i in 1:length(tests_NO)){
      if (is.na(tests_NO[i])){
        tests_NO[i] <- tests_NO[i-1]
      }
    }
    test_t_3 <- tests_NO[(match(as.Date("2020-04-04"), test_timeNO)):(match(as.Date("2021-03-28"), test_timeNO))]
    test_t_4 <- tests_NO[(match(as.Date("2020-04-03"), test_timeNO)):(match(as.Date("2021-03-27"), test_timeNO))]
    test_t_5 <- tests_NO[(match(as.Date("2020-04-02"), test_timeNO)):(match(as.Date("2021-03-26"), test_timeNO))]
    test_t_6 <- tests_NO[(match(as.Date("2020-04-01"), test_timeNO)):(match(as.Date("2021-03-25"), test_timeNO))]
    tests <- cbind(test_t_3, test_t_4, test_t_5, test_t_6)
    
    covariates <- cbind(weighted_mob, tests)
    XX <- t(covariates)%*%covariates
    inv_XX <- chol2inv(chol(XX))
    mean_g <- rep(0, ncol(covariates))
    
    # Modify data into a form suitable for Stan
    reg_data_by_country = list(n_obs = nrow(dependent),
                               K = ncol(dependent),  #number of regressions
                               J = ncol(covariates), #number of covariates
                               x = covariates,
                               y = dependent,
                               b0 = mean_g,
                               inv_XX = inv_XX)
  }
  
  else{
    print("country not included")
  }
  
  return(reg_data_by_country) 
}
