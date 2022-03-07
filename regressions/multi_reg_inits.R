inits_reg <- function(country = c("GR", "PT", "UK", "DE", "SE", "NO")){
  
  if (country == "GR"){ 
    inits_reg_by_country = function(){
      list( delta = rbind(c(runif(1,450,650), runif(1,0.1,0.2),
                             runif(1,-0.025,0.025), runif(1,0.025,0.075), 
                             runif(1, 0.02,0.06)),
                           c(
                             runif(1,0.1410824,0.2691842), runif(1,-0.0000369,-0.0000061), 
                             runif(1,-0.0000252,0.0000180), runif(1,-0.0000284, 0.0000150), 
                             runif(1,-0.0000282,0.0000020)
                           ),
                           c(runif(1,0.1,0.3), runif(1,-0.0001,0.0001), 
                             runif(1,-0.0001,0.0001), runif(1,-0.001,0.001), 
                             runif(1,-0.0001,0.0001))),
             sigma = runif(1,220,265), 
             L_sigma =  c(runif(1,1200,1300), runif(1,1,1.2), runif(1,0.64,0.74))
      )
    }
  }
  
  else if (country == "PT"){
    inits_reg_by_country = function(){
      list( delta = rbind(c(runif(1,450,650), runif(1,0.1,0.2),
                            runif(1,-0.025,0.025), runif(1,0.025,0.075), 
                            runif(1, 0.02,0.06)),
                          c(runif(1,0.1410824,0.2691842), runif(1,-0.0000369,-0.0000061), 
                            runif(1,-0.0000252,0.0000180), runif(1,-0.0000284, 0.0000150), 
                            runif(1,-0.0000282,0.0000020)),
                          c(runif(1,0.1,0.3), runif(1,-0.0001,0.0001), 
                            runif(1,-0.0001,0.0001), runif(1,-0.001,0.001), 
                            runif(1,-0.0001,0.0001))),
             sigma = runif(1,220,265), 
            L_sigma =  c(runif(1,1200,1300), runif(1,1,1.2), runif(1,0.64,0.74))
                    )
                  }
  }
  
  else if (country == "UK"){
    inits_reg_by_country = function(){
      list( delta = rbind(c(runif(1,-11348,-11346), runif(1,0.019,0.021),
                            runif(1,0.01,0.012), runif(1,0.01,0.012), 
                            runif(1, 0.004,0.006)),
                          c(
                            runif(1,0.3,0.4), runif(1,-0.0000005, -0.0000004), 
                            runif(1,-0.00000045,-0.00000035), runif(1,-0.00000045,-0.00000035), 
                            runif(1,-0.0000005, -0.0000004)
                          ),
                          c(runif(1,0.5,0.65), runif(1,-0.00000009,-0.00000007), 
                            runif(1, 0.0000002, 0.0000004), runif(1,0.00000025,0.0000004), 
                            runif(1,0.00000075,0.00000095))),
            sigma = runif(1,400,425), 
            L_sigma =  c(runif(1,2100,2200), runif(1,0.8,1.1), runif(1,0.4,0.5))
      )
    }
  }
  
  else if (country == "DE"){
    inits_reg_by_country = function(){
      list( delta = rbind(c(runif(1,-11350,-11300), runif(1,-0.4,0),
                            runif(1,0.0002,0.0006), runif(1,0.008,0.03), 
                            runif(1, 0.3,0.7)),
                          c(
                            runif(1,0,0.2), runif(1,-0.000012, -0.000008), 
                            runif(1,-0.0000009,-0.0000005), runif(1,-0.0000009,-0.0000005), 
                            runif(1,0,0.000004)
                          ),
                          c(runif(1,0.4,0.6), runif(1,-0.000007,-0.000004), 
                            runif(1,-0.0000006,-0.0000002), runif(1,-0.0000006,-0.0000002), 
                            runif(1,0.000003,0.000007))),
            sigma = runif(1,650,750), 
            L_sigma =  c(runif(1,2300,2800), runif(1,0.6,0.9), runif(1,0.7,1))
      )
    }
  }
  
  else if (country == "SE"){
    inits_reg_by_country = function(){
      list( delta = rbind(c(runif(1,0.1,0.4),runif(1,-0.05,0.05),
                          runif(1,-0.05,0.05), runif(1,-0.05,0.05), 
                          runif(1,-0.05,0.05)),
                        c(runif(1,-0.4,-0.1), runif(1,-0.0001,0.0001),
                          runif(1,-0.0001,0.0001), runif(1,-0.0002,0.0002), 
                          runif(1,-0.0002,0.0002)),
                        c(runif(1,0.1,0.3), runif(1,-0.0001,0.0001), 
                          runif(1,-0.0001,0.0001), runif(1,-0.001,0.001), 
                          runif(1,-0.0001,0.0001))),
           sigma = runif(1,0.5,2), 
           L_sigma =  c(runif(1,700,750), runif(1,0.4,0.7), runif(1,0.4,0.6))
      )
    }
  }
  
  else if (country == "NO"){
    inits_reg_by_country = function(){
      list( delta = rbind(c(runif(1,-55,-45), runif(1,0.0004,0.0008),
                            runif(1,0.001,0.005), runif(1,0.002,0.004), 
                            runif(1, 0.004,0.006)),
                          c(runif(1,-0.1,0), runif(1,-0.00004,-0.00002), 
                            runif(1,-0.00002,0), runif(1,-0.00003, 0), 
                            runif(1,-0.00004,0.00001)),
                          c(runif(1,0.2,0.4), runif(1,0.00003,0.00006), 
                            runif(1,0.00001,0.00004), runif(1,0.00001,0.00004), 
                            runif(1,0.000035,0.00006))),
            sigma = runif(1,30,40), 
            L_sigma =  c(runif(1,95,100), runif(1,0.95,1.1), runif(1,1.1,1.3))
      )
    }
  }
  
  else{
    print("country not included")
  }
  
  return(inits_reg_by_country) 
}