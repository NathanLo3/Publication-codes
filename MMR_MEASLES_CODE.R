# Manuscript: Public health and economic consequences of vaccine hesitancy for measles in the United States
# Nathan C Lo and Peter J Hotez
# JAMA Pediatrics, 2017

# Coded by Nathan C. Lo (Stanford University)
# Email: nathan.lo@stanford.edu 
# Last updated: 4/10/17 

# Model summary 
# This is the R code for the base case analysis in the manuscript. 
# This code requires the associated excel data files.
# 1) CDC_MMR_coverage_state_data.xlsx
# 2) CDC_vacc_exemption_state_data.xlsx
# 3) CDC_MMR_county_data.xlsx
# To reproduce the main analysis, run this script. To change the analysis, skip down to "Run code" section for additional settings!
# **Note: You will need to change your directory files to load these excels.**
setwd("~/Dropbox/Lo collaborations/Peter Hotez/Measles/Model/R/PUBLISHED CODE/") # **Change for your local directory**

# In the main analysis, I simulate ages 2-11 years 
# Primary school population in the United States (kindergarden - 5th grade; ages 5-11). 
# Young children (ages 2-5)

# Load R packages 
library(xlsx)
library(Hmisc)

# Main_or_SA
# 1- Base case (Ages 2-11)
Main_or_SA <- 1 

# Load state-level MMR vaccination data and vaccine exemption data 
data_MMR_coverage <- read.xlsx(file="Data_Model_code_JAMA_Pediatrics_2017.xlsx", sheetName="Data_state", header=TRUE)
data_vacc_exemption <- read.xlsx(file="Data_Model_code_JAMA_Pediatrics_2017.xlsx", sheetName="Data_exemptions", header=TRUE) 

# School VaxView (use 2010-2015, not 2009)
data_MMR_coverage$Coverage_2009 <- data_MMR_coverage$Coverage_2009/100 # convert from percentage (0-100%) to fraction (0-1)
data_MMR_coverage$Coverage_2010 <- data_MMR_coverage$Coverage_2010/100
data_MMR_coverage$Coverage_2011 <- data_MMR_coverage$Coverage_2011/100
data_MMR_coverage$Coverage_2012 <- data_MMR_coverage$Coverage_2012/100
data_MMR_coverage$Coverage_2013 <- data_MMR_coverage$Coverage_2013/100
data_MMR_coverage$Coverage_2014 <- data_MMR_coverage$Coverage_2014/100
data_MMR_coverage$Coverage_2015 <- data_MMR_coverage$Coverage_2015/100

# Child VaxView 
data_MMR_coverage$Coverage_2013b <- data_MMR_coverage$Coverage_2013b/100
data_MMR_coverage$Coverage_2014b <- data_MMR_coverage$Coverage_2014b/100
data_MMR_coverage$Coverage_2015b <- data_MMR_coverage$Coverage_2015b/100

# Load county-data
county_MMR_coverage <- read.xlsx(file="Data_Model_code_JAMA_Pediatrics_2017.xlsx", sheetName="Data_county", header=TRUE)

num_county_perstate <- vector(length=dim(data_MMR_coverage)[1])
SD_perstate <-  vector(length=dim(data_MMR_coverage)[1])
for (i in 1:dim(data_MMR_coverage)[1]) { 
  num_county_perstate[i] <- length(county_MMR_coverage$County[county_MMR_coverage$State==unique(county_MMR_coverage$State)[i]])
  SD_perstate[i] <- sd(county_MMR_coverage$latest_mean[county_MMR_coverage$State==unique(county_MMR_coverage$State)[i]])
}

SD_perstate[is.na(SD_perstate)] <- mean(na.omit(SD_perstate))
SD_perstate <- SD_perstate/100
overall_county_SD <- sd(county_MMR_coverage$latest_mean)/100

# Historical measles outbreak data
# https://www.cdc.gov/measles/cases-outbreaks.html (accessed 3/31/17)
historic_measles_cases <- c(284, 188, 70)
# 2010 2011	2012 2013 2014 2015 2016
# Note 2014 had large outbreak in 2014 in Amish population (~383 cases).
# Overall, from 2013-2015, I estimate approx 180 cases per year

if (Main_or_SA==1) {
  avg_percent_cases <- mean(c((25*0.75 + 36*0.5), (17*0.75 + 25*0.5), (12*0.75 + 17*0.5))) # ages 2-11 (27% of cases)
} 

annual_measles <- mean(historic_measles_cases) * avg_percent_cases/100 
# 2013- 11% <1, 25% 1-4, 36% 5-19
# 2014- 6% <1, 17% 1-4, 25% 5-19
# 2015- 16% <1, 12% 1-4, 17% 5-19 https://www.cdc.gov/mmwr/preview/mmwrhtml/mm6414a1.htm
# Note: I use 50% of 5-19 age bin

annual_3_outbreaks <- mean(c(11,22,15))
num_3_outbreaks <- floor(avg_percent_cases/100*annual_3_outbreaks)
  
# Number of outbreaks >=3
# Average 16 in recent times 

import_cases_per_year <- 34 # 34 import cases/year (min-18, max-80)
# https://www.cdc.gov/vaccines/pubs/surv-manual/chpt07-measles.html (accessed 1/30/17)
# " Of the 911 reported measles cases, 372 (40%) were importations (on average 34 importations/year)"
prob_import_county <- import_cases_per_year/dim(county_MMR_coverage)[1]


# Create US cohort of primary school (ages 5-10; 6 age groups) and pre-school (ages 2-4; 3 age groups) children
# 50 states + DC
# Primary school cohort: Kinder (2015), 1st grade (2014), 2nd grade (2013), 3rd grade (2012), 4th grade (2011), 5th grade (2010)
# Pre-school cohort: 4 years old (2013 child data), 3 years old (2014 child data), 4 years old (2015 child data)
US_children_matrix <- array(dim=c(dim(data_MMR_coverage)[1],4))
for (i in 1:dim(data_MMR_coverage)[1]) { # Loop through 51 states (+DC!)
  
  if (Main_or_SA==1) {
  US_children_matrix[i,1] <- data_MMR_coverage$Coverage_2010[i]*data_MMR_coverage$Kinder_pop_2010[i] + data_MMR_coverage$Coverage_2011[i]*data_MMR_coverage$Kinder_pop_2011[i] + data_MMR_coverage$Coverage_2012[i]*data_MMR_coverage$Kinder_pop_2012[i] + data_MMR_coverage$Coverage_2013[i]*data_MMR_coverage$Kinder_pop_2013[i] + data_MMR_coverage$Coverage_2014[i]*data_MMR_coverage$Kinder_pop_2014[i] + data_MMR_coverage$Coverage_2015[i]*data_MMR_coverage$Kinder_pop_2015[i] + data_MMR_coverage$Coverage_2013b[i]*data_MMR_coverage$Kinder_pop_2015[i] + data_MMR_coverage$Coverage_2014b[i]*data_MMR_coverage$Kinder_pop_2015[i] + data_MMR_coverage$Coverage_2015b[i]*data_MMR_coverage$Kinder_pop_2015[i]
  US_children_matrix[i,2] <- (1-data_MMR_coverage$Coverage_2010[i])*data_MMR_coverage$Kinder_pop_2010[i] + (1-data_MMR_coverage$Coverage_2011[i])*data_MMR_coverage$Kinder_pop_2011[i] + (1-data_MMR_coverage$Coverage_2012[i])*data_MMR_coverage$Kinder_pop_2012[i] + (1-data_MMR_coverage$Coverage_2013[i])*data_MMR_coverage$Kinder_pop_2013[i] + (1-data_MMR_coverage$Coverage_2014[i])*data_MMR_coverage$Kinder_pop_2014[i] + (1-data_MMR_coverage$Coverage_2015[i])*data_MMR_coverage$Kinder_pop_2015[i] + (1-data_MMR_coverage$Coverage_2013b[i])*data_MMR_coverage$Kinder_pop_2015[i] + (1-data_MMR_coverage$Coverage_2014b[i])*data_MMR_coverage$Kinder_pop_2015[i] + (1-data_MMR_coverage$Coverage_2015b[i])*data_MMR_coverage$Kinder_pop_2015[i] 
  } 
}
US_children_matrix[,3] <- (US_children_matrix[,1]+US_children_matrix[,2])
US_children_matrix[,4] <- US_children_matrix[,1]/(US_children_matrix[,1]+US_children_matrix[,2])

US_children_matrix <- data.frame(State=data_MMR_coverage$Names, Kid_pop=US_children_matrix[,3], Coverage=US_children_matrix[,4])
outbreak_risk <- as.integer(US_children_matrix$Coverage < herd_threshold)

US_children_matrix <- data.frame(US_children_matrix, outbreak_risk)
US_children_matrix$State[US_children_matrix$outbreak_risk==1]

county_MMR_coverage$PopSize <- 0
sum_iter3 <- 1
for (j in 1:dim(data_MMR_coverage)[1]) {
  for (k in 1:num_county_perstate[j]) { 
    county_MMR_coverage$PopSize[sum_iter3] <- data_MMR_coverage$Kinder_pop_2015[j]
    sum_iter3 <- sum_iter3+1
  } }

# Main measles model for simulation of outbreaks and cases 
measles_model <- function(sims_model, prob_import_county, R0_mean, hestitent_frac, remove_exemptions) { 

MMR_efficacy <- 0.95 # Vaccine efficacy 

matrix_model <- array(0,dim=c(dim(county_MMR_coverage)[1],sims_model))
for (i in 1:sims_model) {
  sum_iter <- 1
  
    # simulate importation of measles cases 
   for (k in 1:round((dim(county_MMR_coverage)[1]*prob_import_county))) {

      # randomly draw a county for measles introduction, weighted by population size 
      j <- sample(x=seq(1,dim(data_MMR_coverage)[1]), size=1, replace=TRUE, prob=US_children_matrix$Kid_pop/sum(US_children_matrix$Kid_pop))
      
      # simulate MMR coverage in county 
      sim_country_coverage <- rnorm(1, mean=(US_children_matrix$Coverage[j] - hestitent_frac), sd=SD_perstate[j])
      
      # remove exemptions if testing this strategy 
      if (remove_exemptions==1) {
        sim_country_coverage <- sim_country_coverage + data_vacc_exemption$Exemp_2015[j]/100
      }

      # Compute Re for this county 
      R0 <- R0_mean 
      R_eff <- R0* (1-(sim_country_coverage*MMR_efficacy)) # Adds 2008 county-level vaciation to 2015 state estimates
      
      R <- R_eff
      max_outbreak_size <- 100 # Set max outbreak size
      x <- 100
      
      outbreak_store <- vector(length=x)
      outbreak_store_final <- vector(length=x)
      for (q in 1:x) { 
        
        # Stochastic model here (see appendix files for equations)
        num1 <- R^(q-1)
        den1 <- (R+1)^((2*q)-1)
        num2 <- lfactorial((2*q)-2)
        den2 <- (lfactorial(q)+lfactorial(q-1))
        
        outbreak_store[q] <- (num1/den1)*exp(num2-den2)
      }
      
      outbreak_store_final <- cumsum(outbreak_store)
      # 
      
      prob_outbreak_size <- runif(1)
        
        # Compute outbreak size and store
        if (sum(outbreak_store_final >= prob_outbreak_size) > 0) { 
          
          outbreak_size <- min(which(outbreak_store_final >= prob_outbreak_size))
          matrix_model[sum_iter,i] <- outbreak_size
        } else {
          matrix_model[sum_iter,i] <- max_outbreak_size
        }
      #}
      sum_iter <- sum_iter + 1
    }
}

# Summary stats 
mean_num_outbreaks <- mean(colSums(matrix_model>=3))
mean_annual_cases <- mean(colSums(matrix_model))
cat("mean", mean(colSums(matrix_model)), "median", median(colSums(matrix_model)), "lower", quantile(colSums(matrix_model), 0.05), "upper", quantile(colSums(matrix_model), 0.95))
cat("\n", "Number of >=3 outbreaks:", mean(colSums(matrix_model>=3)), "Number of >=1 outbreaks:", mean(colSums(matrix_model>0)))
return(c(median(colSums(matrix_model)), median(colSums(matrix_model>=3)), median(colSums(matrix_model>0)), quantile(colSums(matrix_model), 0.05), quantile(colSums(matrix_model), 0.95)))
}

###############################################################################################
# Measles script: Run code
###############################################################################################
# Here, you can run the code. 
# 1) You can increase prevalence of non-medical exemptions (hestitent_frac_vec). The default
# is to test from 0% increase (current) to 10% increase in each county. 
# 2) You can remove exemptions to childhood vaccination (set "remove_exemptions"=1, default is 0)

sims_model <- 10000 #10,000 runs 
prob_import_county <- 0.06468105 #calibrated 
R0_mean <- 5.67047563 #calibrated
hestitent_frac_vec <- c(seq(0, 0.10, by=0.02)) # You can increase prevalence of non-medical exemptions (hestitent_frac_vec). The default
# is to test from 0% increase (current) to 10% increase in each county. 
cost_per_case <- 20000
remove_exemptions <- 0 # 0- no; 1- yes, remove exemptions 

# measles_model(sims_model, prob_import_county, R0_mean, hestitent_frac, remove_exemptions)

model_output2 <- array(dim=c(length(hestitent_frac_vec), 5))
for (i in 1:length(hestitent_frac_vec)) {

hestitent_frac <- hestitent_frac_vec[i]
model_output2[i,] <- measles_model(sims_model, prob_import_county, R0_mean, hestitent_frac, remove_exemptions)
}

write.csv(model_output2, file="measles_model_output.csv")

####################################################################################
# Plot (data visualization)
####################################################################################

historic_measles_cases2 <- c(63, 220, 55, 187, 284, 188, 70)
hestitent_frac_vec2 <- seq(0,10, by=2)
par(mar=c(5,4,4,5)+.1)
plot(hestitent_frac_vec2, model_output2[,1],col="blue", type="l", yaxs="i",lwd=2, xlim=c(0,10), ylim=c(0,600), ylab="Annual cases", xlab="Non-medical exemptions (%)", cex.axis=1.5, cex.lab=1.25)
polygon( c(hestitent_frac_vec2, rev(hestitent_frac_vec2)), c(model_output2[,4], rev(model_output2[,5])),col=rgb(0,0,1, 0.2),border=NA)
abline(v=2, col="black", type="l", lty=2)
points(rep(2,7), historic_measles_cases2*avg_percent_cases/100 , col=556, bg=556,pch = 21)

par(new=TRUE)
plot(hestitent_frac_vec2, model_output2[,1]*cost_per_case/1000000, pch=., axes=F, col=1, xlab=NA, ylab=NA, cex.axis=1.5, cex.lab=1.25)
axis(side = 4, cex.axis=1.5, cex.lab=1.25)
mtext("Annual costs (million US$)",side = 4, line = 3,cex=1.25)

