# Manuscript: Comparison of strategies and incidence thresholds for Vi conjugate vaccines against typhoid fever: a cost-effectiveness modeling study
# Authored by: Nathan C Lo and colleagues 
# Journal of Infectious Diseases (JID), 2018

# Coded by Nathan C. Lo (Stanford University)
# Email: nathan.lo@stanford.edu 
# Last updated: 11/10/17 

# Model summary 
# This is the R code for the base case analysis in the manuscript. The analysis requires three sequential steps as follows.
# Note:
# 1) Model calibration, use "Run_typhoid_calibration_revised.R" to create posterior distribution of model parameters.
# 2) Model simulation of vaccine strategies, use "Run_typhoid_simulation_revised.R" to create transmission projections
# 3) Cost-effectiveness, use "Run_CEA_typhoid_vaccination.R" to compute cost-effective incidence thresholds

# The CSV age distribution file "Age_distriction_typhoid_111316" is required for models.
# **Note: You will need to change your directory files to load these excels.**

Run_typhoid_simulation <- function(total_incidence) {
  
  file_output_1 <- paste("Transmission_files_1_15_100917_", total_incidence, "_base.csv", sep = '')
  file_output_2 <- paste("Transmission_files_1_15_100917_", total_incidence, "_EPI.csv", sep = '')
  file_output_3 <- paste("Transmission_files_1_15_100917_", total_incidence, "_EPIschool.csv", sep = '')
  file_mcmc_input <- paste("MCMC_chain_091617_", total_incidence, ".csv", sep = '')
  
  Typhoid_model_simulation <- function(total_incidence, file_output_1, file_output_2, file_output_3, file_mcmc_input) { 
    
    ####################################################################################################################################
    #  MODEL: Typhoid age-structured compartmental model with vaccination, SIMULATION template 
    # Date- 
    # Coded by: Nathan C. Lo (Stanford University)
    # Email: Nathan.Lo@stanford.edu 
    
    # Methods:
    # SIRCWRV compartmental model 
    # Time step of 1 month
    # 50 year burn-in period to endemic levels 
    
    # What are we fitting?
    # 1) beta- short cycle transmission coefficient
    # 2) lambda- long cycle transmission coefficient
    # 3) p- reporting/symptomatic fraction
    # 4) w- rate of immune waning from natural immunity to susceptiblity
    # 5) r- relative infectiousness of carrier 
    
    ####################################################################################################################################
    
    # Load R packages  
    require(deSolve) # Differential equation solver 
    #library(xlsx)
    #library(ggplot2)
    library(matrixStats)
    
    # Three setting age-distribution from meta-analysis
    # setwd("~/Documents/R/Typhoid/Transmission_models/ASTMH/") # Use for local machine 
    age_distribution <- read.csv(file="Age_distriction_typhoid_111316.csv", header=FALSE)
    
    # Pick epi #################################################################################
    # total_incidence <- N # Annual incidence: per 100,000
    ############################################################################################
    
    if (total_incidence > 100) {
      age_dist_iter <- age_distribution$V3
    } else if (total_incidence < 10) {
      age_dist_iter <- age_distribution$V1
    } else if (total_incidence>=50 & total_incidence<=100) {
      age_dist_iter <- ((total_incidence-50)/50)*age_distribution$V3 + (1-((total_incidence-50)/50))*age_distribution$V2
    } else if (total_incidence>=10 & total_incidence<50) {
      age_dist_iter <- ((total_incidence-10)/40)*age_distribution$V2 + (1-((total_incidence-10)/40))*age_distribution$V1
    }
    incidence_weight <- age_dist_iter/sum(age_dist_iter)
    incidence_weight <- c(incidence_weight[1]*(1/5), incidence_weight[1]*(4/5), incidence_weight[-1])
    
    # Age structure is 16 age groups (<1, 1-4, 5-9, 10-14, etc until 74)
    age_weight <- c(1.542680, 6.103679, 7.588358, 7.556721, 7.513031, 7.451264, 7.384977, 7.308144, 7.210220, 7.080659, 6.885564, 6.588779, 6.139082, 5.476965, 4.606948, 3.562927)/100
    
    # Morality rate (age-dep), monthly 
    u_age <- c(0.00288262, 0.000188349, 6.68004E-05, 5.00752E-05, 8.35424E-05, 0.000125471, 0.000142272, 0.000167506, 0.000226533, 0.000302733, 0.000439066, 0.000654345, 0.001013536, 0.001663672, 0.002630401, 0.004194527)
    
    # Carriage risk (age-dep)
    carrier_age <-c(0.3, 0.3, 0.3, 0.3, 0.3, 2.1, 2.1, 4.4, 4.4, 8.8, 8.8, 10.1, 10.1, 7.8, 7.8, 7.8)/100
    
    
    total_pop <- 100000*(1/1.25) # Set population at 100,000, with population adjustment to account for demography
    age_group <- length(age_weight)
    # Time step is one month
    
    
    # Compute initial conditions (month, time 0)
    # Keep population fractions (not absolute numbers)
    S.0 <- vector(length=age_group)
    I.0 <- vector(length=age_group)
    for (i in 1:age_group) { 
      I.0[i] <- ((total_incidence*incidence_weight[i])/12)/total_pop # Divide by 12 since converting from annual cases to monthly cases 
      S.0[i] <- (total_pop*age_weight[i] - (I.0[i]*total_pop))/total_pop # If not infected, then they are suceptibles 
    }
    R.0 <- rep(0, age_group)
    C.0 <- rep(0, age_group)
    V.0 <- rep(0, age_group)
    I_in.0 <- rep(0, age_group)
    V_in.0 <- rep(0, 2)
    W.0 <- 0
    
    # Enter initial parameter list 
    # Time step of one month!
    duration_vacc_immune <- 19.2 * 12 # 19.2 years
    
    params <- c(S.0=c(S.0), I.0= c(I.0), R.0= c(R.0), C.0= c(C.0), V.0= c(V.0), W.0=W.0, I_in.0=c(I_in.0), V_in.0=c(V_in.0), lambda=NA, beta=NA,  p=NA, w=NA, r=NA, ksi=NA, w_vacc=1/(duration_vacc_immune), a1=0, a2=0, v=(1/(12*10)), N=total_pop, Bi=(21.2/1000/12), u_i=0.005, gamma=(1/(21/30)), psi=(1/(21/30)))  #  ksi=4.5e9  Vol=(15*30*total_pop)
    # Note, states with multiple transitions (infected -> recovred, carrier) are corrected to ensure average duration is constant 
    
    # Create test typhoid endemic data 
    num_months <- 12
    annual_incidence <- total_incidence
    calibration_years <- 50+10
    incidence1 <- rep(annual_incidence/num_months*incidence_weight[1],num_months*calibration_years)
    incidence2 <- rep(annual_incidence/num_months*incidence_weight[2],num_months*calibration_years)
    incidence3 <- rep(annual_incidence/num_months*incidence_weight[3],num_months*calibration_years)
    incidence4 <- rep(annual_incidence/num_months*incidence_weight[4],num_months*calibration_years)
    incidence5 <- rep(annual_incidence/num_months*incidence_weight[5],num_months*calibration_years)
    incidence6 <- rep(annual_incidence/num_months*incidence_weight[6],num_months*calibration_years)
    incidence7 <- rep(annual_incidence/num_months*incidence_weight[7],num_months*calibration_years)
    incidence8 <- rep(annual_incidence/num_months*incidence_weight[8],num_months*calibration_years)
    incidence9 <- rep(annual_incidence/num_months*incidence_weight[9],num_months*calibration_years)
    incidence10 <- rep(annual_incidence/num_months*incidence_weight[10],num_months*calibration_years)
    incidence11 <- rep(annual_incidence/num_months*incidence_weight[11],num_months*calibration_years)
    incidence12 <- rep(annual_incidence/num_months*incidence_weight[12],num_months*calibration_years)
    incidence13 <- rep(annual_incidence/num_months*incidence_weight[13],num_months*calibration_years)
    incidence14 <- rep(annual_incidence/num_months*incidence_weight[14],num_months*calibration_years)
    incidence15 <- rep(annual_incidence/num_months*incidence_weight[15],num_months*calibration_years)
    incidence16 <- rep(annual_incidence/num_months*incidence_weight[16],num_months*calibration_years)
    
    month <- seq(1,num_months*calibration_years)
    data<- data.frame(month, incidence1, incidence2, incidence3, incidence4, incidence5, incidence6, incidence7, incidence8, incidence9, incidence10, incidence11, incidence12, incidence13, incidence14, incidence15, incidence16)
    
    
    # Compartmental model 
    typhoid.model <- function (t, y, params) { # Note: this input format is required for ode function (used in "prediction" function)
      # x- initial guesses 
      
      age_groups<- age_group # <1, 1-4, 5-9, etc (5 year groups) til 74
      
      S <- vector(length=age_groups)
      I <- vector(length=age_groups)
      R <- vector(length=age_groups)
      C <- vector(length=age_groups)
      V <- vector(length=age_groups)
      W <- vector(length=1)
      
      dS <- vector(length=age_groups)
      dI <- vector(length=age_groups)
      dR <- vector(length=age_groups)
      dC <- vector(length=age_groups)
      dV <- vector(length=age_groups)
      dW <- vector(length=1)
      dI_in <- vector(length=age_groups)
      dV_in <- rep(0, length=2)
      
      for (i in 1:length(y)) {  
        
        if (i<=age_groups) { 
          S[i] <- y[i] 
        } else if (i>age_groups & i<=(2*age_groups)) { 
          I[i-age_groups] <- y[i]
        } else if (i>(2*age_groups) & i<=(3*age_groups)) { 
          R[i-(2*age_groups)] <- y[i]
        } else if (i>(3*age_groups) & i<=(4*age_groups)) { 
          C[i-(3*age_groups)] <- y[i]
        } else if (i>(4*age_groups) & i<=(5*age_groups)) { 
          V[i-(4*age_groups)] <- y[i]
        } else if (i==(5*age_groups + 1)) { 
          W[1] <- y[i]
        } 
      }
      
      N_curr <- (sum(S) + sum(I) + sum(R) + sum(C) + sum(V))*params["N"]
      
      with(
        as.list(params),
        {
          
          for (i in 1:age_groups) { 
            human_contr_temp <- 0
            water_contr_temp <- 0
            for (j in 1:age_groups) { 
              dS[i] <- -beta*S[i]*I[j] - beta*S[i]*C[j]*r + dS[i]
              dI[i] <- beta*S[i]*I[j] + beta*S[i]*C[j]*r + dI[i]
              dR[i] <- dR[i]
              dC[i] <- dC[i]
              dV[i] <- dV[i]
              
              dI_in[i] <- ( beta*S[i]*I[j] + beta*S[i]*C[j]*r)  + dI_in[i] 
              
            }
            
            dS[i] <- -lambda*S[i]*W + w*R[i] + w_vacc*V[i] - u_age[i]*S[i] + dS[i]
            dI[i] <- lambda*S[i]*W  - (1-carrier_age[i])*gamma*I[i] - carrier_age[i]*gamma*I[i] - (u_age[i]+u_i)*I[i] + dI[i]
            dC[i] <- carrier_age[i]*gamma*I[i] - v*C[i] - u_age[i]*C[i] + dC[i]
            dR[i] <- (1-carrier_age[i])*gamma*I[i] - w*R[i] + v*C[i] - u_age[i]*R[i] + dR[i]
            dV[i] <- -w_vacc*V[i] - u_age[i]*V[i] + dV[i]
            # dW <- ksi*I[i]  + ksi*r*C[i] - psi*W + dW
            dW <- ksi*I[i]  + ksi*r*C[i] + dW
            if (i==1) { 
              dW <- -psi*W + dW
            }
            
            dI_in[i]  <- lambda*S[i]*W + dI_in[i] 

            # Enter aging process
            # Age structure (16 groups)
            # <1 yr
            # 1-4 year
            # 5-9 year
            # 10-14 yr
            # etc
            # max age, 74 yrs
            
            if (i==1) { 
              # Age group <1
              # 1) Model births. Add to suceptibles compartment at birth rate.
              # 2) Age up the infant cohort; 1/12 per month since 1 year group
              dS[i] <- dS[i] + Bi - (1/12)*S[i]
              dI[i] <- dI[i] - (1/12)*I[i]
              dR[i] <- dR[i] - (1/12)*R[i]
              dC[i] <- dC[i] - (1/12)*C[i]
              dV[i] <- dV[i] - (1/12)*V[i]
              
            } else if (i==2) { 
              # Age group 1-4
              # 1) Age upthe cohort (per year) in all compartments; 1/48 per month since 4 year group
              dS[i] <- dS[i] + (1/12)*S[i-1] - (1/48)*S[i]
              dI[i] <- dI[i] + (1/12)*I[i-1] - (1/48)*I[i]
              dR[i] <- dR[i] + (1/12)*R[i-1] - (1/48)*R[i]
              dC[i] <- dC[i] + (1/12)*C[i-1] - (1/48)*C[i]
              dV[i] <- dV[i] + (1/12)*V[i-1] - (1/48)*V[i]
              
            } else if (i==3) { 
              # Age group 5-9
              # 1) Age upthe cohort (per year) in all compartments; 1/48 per month since 4 year group
              dS[i] <- dS[i] + (1/48)*S[i-1] - (1/60)*S[i]
              dI[i] <- dI[i] + (1/48)*I[i-1] - (1/60)*I[i]
              dR[i] <- dR[i] + (1/48)*R[i-1] - (1/60)*R[i]
              dC[i] <- dC[i] + (1/48)*C[i-1] - (1/60)*C[i]
              dV[i] <- dV[i] + (1/48)*V[i-1] - (1/60)*V[i]
              
            } else  { 
              # Age groups (5-year categories) from 10 until 74. Then all die (yikes)
              # 1) Age upthe cohort (per year) in all compartments; 1/60 per month since 5 year group
              dS[i] <- dS[i] + (1/60)*S[i-1] - (1/60)*S[i]
              dI[i] <- dI[i] + (1/60)*I[i-1] - (1/60)*I[i]
              dR[i] <- dR[i] + (1/60)*R[i-1] - (1/60)*R[i]
              dC[i] <- dC[i] + (1/60)*C[i-1] - (1/60)*C[i]
              dV[i] <- dV[i] + (1/60)*V[i-1] - (1/60)*V[i]
            }
            
            # EPI vaccination: ages <1 group only 
            # Note: Infected, carriers, and water (Directly) are not relevant here.
            if (i==1) { # For age group <1 
              if (a1>0) { # when we vaccinate, this will be non-zero
                dS[i] <- -Bi*a1 + dS[i] # vaccinate *fraction* of newborns (remove from suceptible, since introduced in aging section)
                dV[i] <- Bi*a1 + dV[i] # vaccinate *fraction* of newborns
                
                dV_in[1] <- Bi*a1 + dV_in[1]
              } 
            }
            
            # Vaccination of school-aged children: ages 5-9 and 10-14
            # Note: Infected, carriers, and water (Directly) are not relevant here.
            # One month pulse (convert to rate)
            
            # Ages 1-15 years i=2, 3, 4
            # Ages 1-30 years i=2, 3, 4, 5, 6, 7
            
            if (i>=2 & i<=4) {
              
              a2_rate <- -log(1-a2)/1
              dS[i] <- -a2_rate*S[i] + dS[i] # Vaccinate suceptibles and recovered at equal rate
              dV[i] <- a2_rate*S[i] + dV[i] # + a2_calc*R[i]
              
              dV_in[2] <- a2_rate*S[i] + a2*C[i] + a2*R[i]  + dV_in[2] # Note rate for dynamic susceptibles, fraction for stables (recovered, carrier)
            }
            
          }
          dx<- c(dS, dI, dR, dC, dV, dW, dI_in, dV_in)
          list(dx)
          #print(t)
        }
      )
    }
    
    prediction <- function (params, times) {
      xstart <- params[c("S.01","S.02","S.03","S.04","S.05","S.06","S.07","S.08","S.09","S.010","S.011","S.012","S.013","S.014","S.015","S.016","I.01","I.02","I.03","I.04","I.05","I.06","I.07","I.08","I.09","I.010","I.011","I.012","I.013","I.014","I.015", "I.016","R.01","R.02","R.03","R.04","R.05","R.06","R.07","R.08","R.09","R.010","R.011","R.012","R.013","R.014","R.015","R.016","C.01","C.02","C.03","C.04","C.05","C.06","C.07","C.08","C.09","C.010","C.011","C.012","C.013","C.014","C.015","C.016","V.01","V.02","V.03","V.04","V.05","V.06","V.07","V.08","V.09","V.010","V.011","V.012","V.013","V.014","V.015","V.016","W.0", "I_in.01", "I_in.02", "I_in.03", "I_in.04", "I_in.05", "I_in.06", "I_in.07", "I_in.08", "I_in.09", "I_in.010", "I_in.011", "I_in.012", "I_in.013", "I_in.014", "I_in.015", "I_in.016", "V_in.01", "V_in.02")]  # ODE function 
      # y- initial conditions
      # func- set of ODEs
      # params- passed to function 
      out <- ode(
        func=typhoid.model,
        y=xstart,
        times= c(0, times),
        parms=params,
        method="rk4"
      )
      tail(out[, 83:100], num_months*10) # return the I variable only
      #out <- data.frame(out)
    }
    
    prediction2 <- function (params, times) {
      xstart <- params[c("S.01","S.02","S.03","S.04","S.05","S.06","S.07","S.08","S.09","S.010","S.011","S.012","S.013","S.014","S.015","S.016","I.01","I.02","I.03","I.04","I.05","I.06","I.07","I.08","I.09","I.010","I.011","I.012","I.013","I.014","I.015", "I.016","R.01","R.02","R.03","R.04","R.05","R.06","R.07","R.08","R.09","R.010","R.011","R.012","R.013","R.014","R.015","R.016","C.01","C.02","C.03","C.04","C.05","C.06","C.07","C.08","C.09","C.010","C.011","C.012","C.013","C.014","C.015","C.016","V.01","V.02","V.03","V.04","V.05","V.06","V.07","V.08","V.09","V.010","V.011","V.012","V.013","V.014","V.015","V.016","W.0", "I_in.01", "I_in.02", "I_in.03", "I_in.04", "I_in.05", "I_in.06", "I_in.07", "I_in.08", "I_in.09", "I_in.010", "I_in.011", "I_in.012", "I_in.013", "I_in.014", "I_in.015", "I_in.016", "V_in.01", "V_in.02")]  # ODE function 
      # func- set of ODEs
      # params- passed to function 
      out <- ode(
        func=typhoid.model,
        y=xstart,
        times= c(0, times),
        parms=params,
        method="rk4" 
      )
      out <- data.frame(out)
    }
    
    # Define likelihood function with poisson distribution (count data)
    poisson.loglik <- function (params, data) {
      times <- data$month # convert to time-scale year
      pred <- prediction(params,times)*params["p"]*params["N"]
      
      pred_I2 <- matrix(nrow=dim(pred)[1], ncol=dim(pred)[2])
      pred_I2[1,] <- rep(0,dim(pred)[2])
      for (j in 2:dim(pred)[1]) {
        pred_I2[j,] <- as.numeric(pred[j,] - pred[j-1,])
      }
      
      pred_I <- matrix(nrow=(dim(pred)[1]/12), ncol=dim(pred)[2])
      for (j in 1:(dim(pred)[1]/12)) {
        pred_I[j,] <- colSums(pred_I2[((j-1)*12+1):(j*12),])
      }
      pred_I[1,] <- pred_I[1,]/(11/12)
      
      # For MCMC
      # Use subset of age groups (use >1 case per month)
      x1 <- sum(dpois(x=(round(tail(data$incidence1*12, 10))+1),lambda=((pred_I[,1])+1), log=TRUE))
      x2 <- sum(dpois(x=(round(tail(data$incidence2*12, 10))+1),lambda=((pred_I[,2])+1), log=TRUE))
      x3 <- sum(dpois(x=(round(tail(data$incidence3*12, 10))+1),lambda=((pred_I[,3])+1), log=TRUE))
      x4 <- sum(dpois(x=(round(tail(data$incidence4*12, 10))+1),lambda=((pred_I[,4])+1), log=TRUE))
      x5 <- sum(dpois(x=(round(tail(data$incidence5*12, 10))+1),lambda=((pred_I[,5])+1), log=TRUE))
      x6 <- sum(dpois(x=(round(tail(data$incidence6*12, 10))+1),lambda=((pred_I[,6])+1), log=TRUE))
      x7 <- sum(dpois(x=(round(tail(data$incidence7*12, 10))+1),lambda=((pred_I[,7])+1), log=TRUE))
      x8 <- sum(dpois(x=(round(tail(data$incidence8*12, 10))+1),lambda=((pred_I[,8])+1), log=TRUE))
      x9 <- sum(dpois(x=(round(tail(data$incidence9*12, 10))+1),lambda=((pred_I[,9])+1), log=TRUE))
      x10<- sum(dpois(x=(round(tail(data$incidence10*12, 10))+1),lambda=((pred_I[,10])+1), log=TRUE))
      x11<- sum(dpois(x=(round(tail(data$incidence11*12, 10))+1),lambda=((pred_I[,11])+1), log=TRUE))
      x12<- sum(dpois(x=(round(tail(data$incidence12*12, 10))+1),lambda=((pred_I[,12])+1), log=TRUE))
      x13<- sum(dpois(x=(round(tail(data$incidence13*12, 10))+1),lambda=((pred_I[,13])+1), log=TRUE))
      x14<- sum(dpois(x=(round(tail(data$incidence14*12, 10))+1),lambda=((pred_I[,14])+1), log=TRUE))
      x15<- sum(dpois(x=(round(tail(data$incidence15*12, 10))+1),lambda=((pred_I[,15])+1), log=TRUE))
      x16<- sum(dpois(x=(round(tail(data$incidence16*12, 10))+1),lambda=((pred_I[,16])+1), log=TRUE))

      totes <- sum(x1 , x2 , x3 , x4 , x5 , x6 , x7 , x8 , x9 , x10 , x11 , x12 , x13 , x14 , x15, x16)
      return(totes)
    }
    
    # f
    # Objective: How does a set of model parameters perform for poisson negative log-likelihood?
    likelihood <- function (par_initial) {
      
      par <- params
      par[c("beta")] <- par_initial[1]
      par[c("p")] <- par_initial[2]
      par[c("lambda")] <- par_initial[3]
      par[c("w")] <- par_initial[4]
      par[c("r")] <- par_initial[5]
      par[c("ksi")] <- par_initial[6]
      cat("beta",par[c("beta")], "lambda", par[c("lambda")], "p",par[c("p")], "w",par[c("w")], "r",par[c("r")], "ksi",par[c("ksi")])
      poiss <- poisson.loglik(par,data) # smaller numer is better fit (20 better than 200); should find global min
      print(poiss)
      poiss # min with negative 
    }
    
    
    
    # ####################################################################################################################################
    # 
    # # Simulated interventions 
    # ####################################################################################################################################
    # setwd("~/Dropbox/Lo Andrews Shared/Ongoing projects/Typhoid model/Code store/Calibration/")  # Use for local machine 
    
    mcmc_parameters <- read.csv(file=file_mcmc_input, header=TRUE)
    mcmc_parameters <- mcmc_parameters[,-1]
    colnames(mcmc_parameters)<- c("beta", "p","lambda", "w", "r", "ksi", "likelihood")
    burnIn = 1000 # already burned in
    posterior_accepted <- mcmc_parameters[-(1:burnIn),]
    posterior_accepted$prob <- exp(posterior_accepted$likelihood)
    
    # Base case estimate
    #base_case_parameters <- posterior_accepted[posterior_accepted$likelihood==max(posterior_accepted$likelihood),-6][1,]
    vacc_coverage_EPI <- 0.85 
    vacc_coverage_school <- 0.75 
    vacc_efficacy <- 0.915
    
    total_time <- 10*12
    
    simulation_size <- 1000
    transmission_store_matrix1 <- array(0, c(total_time, 18, simulation_size))
    transmission_store_matrix2 <- array(0, c(total_time, 18, simulation_size))
    transmission_store_matrix3 <- array(0, c(total_time, 18, simulation_size))
    
    for (k in 1:simulation_size) {
      # PSA
      # iter_parameters <- sample(x=seq(1,dim(posterior_accepted)[1]), size=1, replace=TRUE, prob=posterior_accepted$prob)
      iter_parameters <- sample(x=seq(1,dim(posterior_accepted)[1]), size=1, replace=TRUE) # no weight, treat all equally
      
      params["beta"] <- posterior_accepted[iter_parameters,1]
      params["lambda"] <- posterior_accepted[iter_parameters,3]
      params["p"] <- posterior_accepted[iter_parameters,2]
      params["w"] <- posterior_accepted[iter_parameters,4]
      params["r"] <- posterior_accepted[iter_parameters,5]
      params["ksi"] <- posterior_accepted[iter_parameters,6]
      
      
      # Run model with fitted parameters and burn-in for 50 yeras
      
      # Burn-in period 50 years, obtain inputs
      time_burn <- 50
      times <- seq(1,(time_burn*12))
      out_calibration<- data.frame(prediction2(params,times))
      initial_calibrated <- tail(out_calibration[2:(length(out_calibration))],1) # Pull new initial conditions from quasi-equil state @ time=50 years
      params_calibrated <- params
      params_calibrated[1:81] <- as.numeric(initial_calibrated[1:81])
      
      ####################################################################################################################################
      # Tested vaccination strategies with calibrated model!
      # Strategies (10 year simulation)
      # 1) model.pred_basecase - base case, no intervention
      # 2) EPI + school catch-up campaign - 85% coverage, instantaneous, <1 and 5-14 year olds
      # 3) EPI - 85% coverage, instantaneous, <1  year olds
      #total_time <- 12*10
      
      # 1) Base case analysis (model.pred_basecase)
      times_intervention <- seq(1, total_time) # X year simulation (one month time step)
      model.pred_basecase <- data.frame(round(prediction2(params_calibrated,times_intervention)[,83:100]*params["N"])) * params["p"]
      
      # base case
      pred_I_vacc_base <- matrix(nrow=dim(model.pred_basecase)[1], ncol=dim(model.pred_basecase)[2])
      pred_I_vacc_base[1,] <- rep(0,dim(model.pred_basecase)[2])
      for (j in 2:dim(model.pred_basecase)[1]) {
        pred_I_vacc_base[j,] <- as.numeric(model.pred_basecase[j,] - model.pred_basecase[j-1,])
      }
      pred_I_vacc_base <- pred_I_vacc_base[-1,]
      
      p_old <- params["p"]
      params["p"] <- (total_incidence/(sum(pred_I_vacc_base)/10))*params["p"]
      pred_I_vacc_base <- pred_I_vacc_base*(params["p"]/p_old)
      
      ###############################################################################################################################
      # Intervention
      # a1- EPI
      # a2- school catch-up
      
      # 1) Vaccine intervention #1 (EPI + school catch-up campaign)
      
      # One-time mass catch-up vaccination of school-aged children (ages 5-14) with EPI
      time_intervention <- 1 # Catch-up school campaign in population of 100,000 done in one month
      #vacc_coverage_EPI <- 0.85 # EPI coverage estimation
      #vacc_coverage_school <- 0.75 # EPI coverage estimation
      #vacc_efficacy <- 0.915
      params_calibrated["a1"] <- vacc_coverage_EPI*vacc_efficacy
      params_calibrated["a2"] <- vacc_coverage_school*vacc_efficacy 
      
      
      times_vaccine <- seq(1,time_intervention) #12 months
      model.pred_vaccine1 <- data.frame(prediction2(params_calibrated,times_vaccine))
      
      # Model EPI, beyond the one-time school catch up
      time_post_vaccine <- seq(1,((total_time)-time_intervention)) # Total of 10 year (minus 1 month) simulation
      params_calibrated["a1"] <- vacc_coverage_EPI*vacc_efficacy #
      params_calibrated["a2"] <- 0 # Turn off school-based vaccine (no intervention)
      params_post_vaccine1 <- tail(model.pred_vaccine1[2:(length(model.pred_vaccine1))],1) # Pull new initial conditions from quasi-equil state post-intervention
      params_calibrated2 <- params_calibrated
      params_calibrated2[1:99] <- as.numeric(params_post_vaccine1[1:99])
      model.pred_vaccine1_part2 <- data.frame(prediction2(params_calibrated2,time_post_vaccine)[,83:100]*params["N"])
      model_vaccine1_full <- rbind((model.pred_vaccine1[,83:100]*params["N"]), model.pred_vaccine1_part2)
      model_vaccine1_full <- rbind(model_vaccine1_full[1:dim(model.pred_vaccine1)[1],], model_vaccine1_full[(dim(model.pred_vaccine1)[1]+2):dim(model_vaccine1_full)[1],])
      model_vaccine1_full <- model_vaccine1_full *params["p"]
      
      model_vaccine1_full[,c(17,18)] <- model_vaccine1_full[,c(17,18)]/params["p"]
      
      #plot(model_vaccine1_full$I.01)
      
      # vaccine
      pred_I_vacc_SAC_EPI <- matrix(nrow=dim(model_vaccine1_full)[1], ncol=dim(model_vaccine1_full)[2])
      pred_I_vacc_SAC_EPI[1,] <- rep(0,dim(model_vaccine1_full)[2])
      for (j in 2:dim(model_vaccine1_full)[1]) {
        pred_I_vacc_SAC_EPI[j,] <- as.numeric(model_vaccine1_full[j,] - model_vaccine1_full[j-1,])
      }
      pred_I_vacc_SAC_EPI <- pred_I_vacc_SAC_EPI[-1,]
      
      
      ###############################################################################################################################
      # Intervention
      
      # 2) Vaccine intervention #2 (EPI only)
      # Keeling Rohani 8.1.1 Pediatric Vaccination
      
      # Annual EPI vaccination of ages 0-4
      time_intervention <- total_time
      #vacc_coverage_EPI <- 0.85 # EPI coverage estimation
      #vacc_coverage_school <- 0.75 # EPI coverage estimation
      #vacc_efficacy <- 0.915
      params_calibrated["a1"] <-  vacc_coverage_EPI*vacc_efficacy 
      params_calibrated["a2"] <- 0
      
      times_vaccine <- seq(1,time_intervention) #1 month
      model.pred_EPI <- data.frame(prediction2(params_calibrated,times_vaccine)[,83:100]*params["N"]*params["p"])
      model.pred_EPI[,c(17,18)] <- model.pred_EPI[,c(17,18)]/params["p"]
      
      # vaccine
      pred_I_vacc_EPI <- matrix(nrow=dim(model.pred_EPI)[1], ncol=dim(model.pred_EPI)[2])
      pred_I_vacc_EPI[1,] <- rep(0,dim(model.pred_EPI)[2])
      for (j in 2:dim(model.pred_EPI)[1]) {
        pred_I_vacc_EPI[j,] <- as.numeric(model.pred_EPI[j,] - model.pred_EPI[j-1,])
      }
      pred_I_vacc_EPI <- pred_I_vacc_EPI[-1,]
      
      
      # Plot the figure
      # plot(seq(0,total_time/12 - 1/12, 1/12), rep(mean(rowSums(pred_I_vacc_base)),dim(pred_I_vacc_base)[1]) , col="blue", type="l", lwd=2, ylim=c(0,50), ylab="Monthly incidence", xlab="Year", cex.axis=1.5, cex.lab=1.5)
      # lines(seq(0,total_time/12 - 1/12, 1/12), rowSums(pred_I_vacc_SAC_EPI), type="l", lwd=2,  col="palegreen4")
      # lines(seq(0,total_time/12 - 1/12, 1/12), rowSums(pred_I_vacc_EPI), type="l", lwd=2,  col="red")
      # legend(20, 18, c("Base case", "EPI", "EPI+School"), xjust = 0.5, lty=c(1,1), lwd=c(2.5,2.5), xpd = TRUE, col=c("blue", "red","palegreen4"), inset = c(0,0), bty="n")
      
      transmission_store_matrix1[,,k] <- pred_I_vacc_base
      transmission_store_matrix2[,,k] <- pred_I_vacc_EPI
      transmission_store_matrix3[,,k] <- pred_I_vacc_SAC_EPI
    }
    
    transmission_summary1 <- array(0, c(total_time,17*2))
    transmission_summary2 <- array(0, c(total_time,17*2))
    transmission_summary3 <- array(0, c(total_time,17*2))
    
    for (i in 1:18) {
      
      if (i <= 16) { 
        transmission_summary1[,1+2*(i-1)] <- rowMeans(transmission_store_matrix1[,i,])
        transmission_summary1[,(i*2)] <- rowSds(transmission_store_matrix1[,i,])
        
        transmission_summary2[,1+2*(i-1)] <- rowMeans(transmission_store_matrix2[,i,])
        transmission_summary2[,(i*2)] <- rowSds(transmission_store_matrix2[,i,])
        
        transmission_summary3[,1+2*(i-1)] <- rowMeans(transmission_store_matrix3[,i,])
        transmission_summary3[,(i*2)] <- rowSds(transmission_store_matrix3[,i,])
      } else {
        transmission_summary1[,(16+i)] <- rowMeans(transmission_store_matrix1[,i,])
        transmission_summary2[,(16+i)] <- rowMeans(transmission_store_matrix2[,i,])
        transmission_summary3[,(16+i)] <- rowMeans(transmission_store_matrix3[,i,])
      }
    }
    
    
    
    write.csv(transmission_summary1, file= file_output_1)
    write.csv(transmission_summary2, file= file_output_2)
    write.csv(transmission_summary3, file= file_output_3)
    
  }
  # Run the typhoid model 
  Typhoid_model_simulation(total_incidence, file_output_1, file_output_2, file_output_3, file_mcmc_input) 
  
}