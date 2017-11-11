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

library(xlsx)
setwd("~/Dropbox/Lo Andrews Shared/Ongoing projects/Typhoid model/Code store/FINAL/Template_final/revise/")

ICER_compute <- function(incidence_scenario, mortality_rate, sim_years, WTPT) { 
  
  Societal=TRUE
  
  # Insert transmission files name here.
  file_base <- paste("Transmission_files_092717_", incidence_scenario, "_base.csv", sep = '')
  file_EPI <- paste("Transmission_files_092717_", incidence_scenario, "_EPI.csv", sep = '')
  file_EPIschool <- paste("Transmission_files_092717_", incidence_scenario, "_EPIschool.csv", sep = '')
  
  data_base=read.csv(file_base)[1:(sim_years*12),-1]
  data_EPI=read.csv(file_EPI)[1:(sim_years*12),-1]
  data_school_EPI=read.csv(file_EPIschool)[1:(sim_years*12),-1]

  # Who was vaccinated
  data_EPI_vaccinated <- data_EPI[,c(33,34)]
  data_school_EPI_vaccinated <- data_school_EPI[,c(33,34)]
  
  # Store SD
  data_base_SD <- data_base[,seq(2,32,by=2)]
  data_EPI_SD <- data_EPI[,seq(2,32,by=2)]
  data_school_EPI_SD <- data_school_EPI[,seq(2,32,by=2)]
  
  # Keep prevalence alone 
  data_base <- data_base[,seq(1,31,by=2)]
  data_EPI <- data_EPI[,seq(1,31,by=2)]
  data_school_EPI <- data_school_EPI[,seq(1,31,by=2)]

age_groups <- 16

# Cost 
vacc_cost <- 2.50 # 2.50 per dose (one dose) gavi price 

EPI_per_person_cost <- (1.10) + vacc_cost  # EPI (1.10)
school_per_person_cost <- (1.30) + vacc_cost # School-based vaccination (1.30)

cost_EPI_annual <- EPI_per_person_cost # age_weight[1] is age 0-4
cost_school_annual <- school_per_person_cost

cost_matrix_base <- rep(0,sim_years)

cost_matrix_EPI <- vector(length=sim_years)
cost_matrix_school_EPI <- vector(length=sim_years)
for (k in 1:sim_years) {
  cost_matrix_EPI[k] <- sum(data_EPI_vaccinated[,1][((k-1)*12+1):(12*k)])*EPI_per_person_cost
  cost_matrix_school_EPI[k] <- sum(data_school_EPI_vaccinated[,1][((k-1)*12+1):(12*k)])*EPI_per_person_cost + sum(data_school_EPI_vaccinated[,2][((k-1)*12+1):(12*k)])*school_per_person_cost
}

cost_matrix <- data.frame(cost_matrix_base, cost_matrix_EPI, cost_matrix_school_EPI)

child_cost_illness <- (0.09*222+0.91*33)
adult_cost_illness <- (0.27*444+0.73*66)

societal_cost <- array(0, dim=c(sim_years, 3))
for (j in 1:sim_years) {
  low <-  1 + 12*(j-1)
  high <- 12*j
  
  societal_cost[j,1] <- sum((data_base[low:high,1:4])*child_cost_illness) + sum((data_base[low:high,5:16])*adult_cost_illness)
  societal_cost[j,2] <- sum((data_EPI[low:high,1:4])*child_cost_illness) + sum((data_EPI[low:high,5:16])*adult_cost_illness)
  societal_cost[j,3] <- sum((data_school_EPI[low:high,1:4])*child_cost_illness) + sum((data_school_EPI[low:high,5:16])*adult_cost_illness)
  
}

if (Societal==TRUE) {
  cost_matrix[,1:3] <- cost_matrix[,1:3] + societal_cost
}

# Disability 

DALY_matrix <- array(dim=c(120,6,3)) # 120 months, 6 disability (4 states, one death, total), and 3 scenarios (base case, EPI, EPI+school)

# Mortality (by age and discounted)
YLL_age_stratified <- c(69, 70.5, 67.1, 62.4, 57.55, 52.9, 48.25, 43.65, 39.05, 34.55, 30.15, 25.9, 21.8, 18, 14.55, 11.55) # WHO Life Tables (South Asia)
# http://apps.who.int/gho/data/view.main.LIFESEAR?lang=en accessed 1/31/17
disc <- 0.03

DALY_mortality <- array(0, dim=c(age_groups, sim_years)) # 16 age groups by 10 year simulation 
for (i in 1:sim_years) { # 16 age groups 
  
  for (j in 1:age_groups) { # 10 year simulation 
    LE <- round(YLL_age_stratified[j])
    
    for (k in i:(LE+i-1)) { # discount based on LE at age 
      discount <- 1/(1+disc)^(k-1)
      DALY_mortality[j,i] <- DALY_mortality[j,i] + discount
    }
  }
}

# Disability 
#Disability Weight Associated with each state
wGI <- 0.325
wIP <- 0.324
wSevere <- 0.133
wMild <- 0.006

#Proportions
rGI <- .17
rIP <- .0025
rSevere <- .4775
rMild <- .35

# Time
duration <- 3/52 # 3 weeks duration of infection 

# Compute disability 
base_disability_matrix <- rowSums(data_base*duration* (rGI*wGI + rIP*wIP + rSevere*wSevere + rMild*wMild))
EPI_disability_matrix <- rowSums(data_EPI*duration* (rGI*wGI + rIP*wIP + rSevere*wSevere + rMild*wMild))
school_EPI_disability_matrix <- rowSums(data_school_EPI*duration* (rGI*wGI + rIP*wIP + rSevere*wSevere + rMild*wMild))

base_disability <- vector(length=sim_years)
EPI_disability <- vector(length=sim_years)
school_EPI_disability <- vector(length=sim_years)

for (i in 1:sim_years) {
  base_disability[i] <- sum(base_disability_matrix[((i-1)*12+1):((i)*12)])
  EPI_disability[i] <- sum(EPI_disability_matrix[((i-1)*12+1):((i)*12)])
  school_EPI_disability[i] <- sum(school_EPI_disability_matrix[((i-1)*12+1):((i)*12)])
}



# Compute mortality 
mortality_matrix <- matrix(0, nrow=3, ncol=age_groups)
for (i in 1:age_groups) {
  for (j in 1:sim_years) { 
    low <-  1 + 12*(j-1)
    high <- 12*j
    
  mortality_matrix[1,i] <- mortality_matrix[1,i] + sum(data_base[low:high,i] * mortality_rate * DALY_mortality[i,j])
  mortality_matrix[2,i] <- mortality_matrix[2,i] + sum(data_EPI[low:high,i] * mortality_rate * DALY_mortality[i,j])
  mortality_matrix[3,i] <- mortality_matrix[3,i] + sum(data_school_EPI[low:high,i] * mortality_rate * DALY_mortality[i,j])
  }
}
mortality_sum <- rowSums(mortality_matrix)


# Discount cost and disability
# Note: mortality is already discounted, so don't include here. 
  for (i in 1:sim_years) {
    
    discount <- 1/(1+disc)^(i-1)
    
    # Disability
    base_disability[i] <-  base_disability[i] * discount
    EPI_disability[i] <- EPI_disability[i] * discount
    school_EPI_disability[i] <- school_EPI_disability[i] * discount
    
    # Cost
    cost_matrix[i,] <- cost_matrix[i,] * discount
  }
  

# Compute ICER

total_cost <- as.vector(colSums(cost_matrix))
total_disability <- c((sum(base_disability)+mortality_sum[1]), (sum(EPI_disability)+mortality_sum[2]), (sum(school_EPI_disability)+mortality_sum[3]))
  

ICER1 <- (total_cost[2] - total_cost[1]) / (total_disability[1] - total_disability[2])
ICER2 <- (total_cost[3] - total_cost[2]) / (total_disability[2] - total_disability[3])

ICER_output <- c(ICER1, ICER2)
return(ICER_output)
} 

mortality_rate_vec <- c(0.005) # Base case of 0.5% CFR

# WTPT
WTPT <- 1035
  

ICER_tracker <- array(dim=c(length(mortality_rate_vec),2))
for (k in 1:length(mortality_rate_vec)) { 
  
mortality_rate <- mortality_rate_vec[k] # 0.01 is 1% mortality 

sim_years <- 10 # 120 iterations (12 month per year)
incidence_scenario_vector <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000)
ICER_store <- array(0, dim=c(length(incidence_scenario_vector),2))
for (i in 1:length(incidence_scenario_vector)) { 
incidence_scenario <- incidence_scenario_vector[i]
ICER_output <- ICER_compute(incidence_scenario, mortality_rate, sim_years, WTPT)

ICER_store[i,] <- ICER_output
}

incidence_scenario_vector_iterp <- incidence_scenario_vector

# linear interp
threshold_EPI <- approx(y=incidence_scenario_vector_iterp, x=ICER_store[,1], xout=WTPT, method="linear")$y
threshold_EPI_school <- approx(y=incidence_scenario_vector_iterp, x=ICER_store[,2], xout=WTPT, method="linear")$y

cat("Threshold for EPI:", threshold_EPI, "Threshold for EPI+school:", threshold_EPI_school)

ICER_tracker[k,] <- c(threshold_EPI, threshold_EPI_school)

}



