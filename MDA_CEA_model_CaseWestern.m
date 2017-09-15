function [total_DALY_sequelae_pre, total_DALY_sequelae_c, total_DALY_sequelae_a, ICER_values, total_DALYs_age, total_cost, total_DALY_iteration] = MDA_CEA_model_CaseWestern(MDA_strategy, snail_strategy, g_MDA, g_snail, MDA_freq, snail_freq, excel_decision)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cost-effectiveness analysis for control/elimination of schistosomiasis
% We test MDA, molluscicides, integrated MDA+molluscicides
% Adapted from Lo et al (2016) Lancet Infect Dis
% Written by Nathan Lo
%
% Last updated: 8/16/17
%
% Collaboration with Case Western
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY
% This function computes cost-effectiveness of schistosomiasis
% control and elimination strategies. We use results from a dynamic
% transmission model (human-snail coupled model, Case Western) that
% simulates various interventions. We input prevalence of no, light, and
% heavy infections for S. haematobium for each age sub-group (pre-school
% aged children, school-aged children, and adults). We then compute the
% total costs (2016 USD), disability (DALYs), and the ICER values for
% specified straegies including MDA,  molluscicides, and an integrated
% MDA+molluscicide strategy.This simulation is a a 10-year period.
% Disability is assigned following prior work and as outlined in the
% Methods/appendix. Disability is reported by age, sequelae, and total.

% The transmission modeling is provided by David Gurarie and Charles King
% (Case Western), and results are shared through excel files.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%
% 1) MDA_strategy (scalar)
%        (0)None
%        (1)MDA (SAC)
%        (2)MDA (SAC, pre-SAC)
%        (3)MDA (Community)
% 2) Snail_strategy (scalar)
%        (0)None
%        (1)Molluscicides
% 3) MDA_coverage (g_MDA; 1x3 vector; fraction e.g., 75% is "0.75")
%        Note: [pre-SAC, SAC, adults]
% 4) Snail_coverage (g_snail; scalar; fraction e.g., 75% is "0.75") )
% 5) MDA frequency (MDA_freq)
%        (1)Annual
%        (2)Biannual (2/year)
% 6) Molluscicides frequency (snail_freq)
%        (1)Annual
%        (2)Biannual (2/year)
%        (12)Monthly
% 7) excel_decision
%        (1) Which sheet (1, 2, 3...)
%        (2) Base case (0), lower limit (-1), upper limit (1)

% Each input matrix is 12x121. The columns represent time (121). Each
% column is the prevalence measurement at baseline (time_0) and for 10
% years at 1 month time steps. The rows represent age-specific prevalence
% of light and heavy infection (mean and SD) from simulation.
% Note: Focus on five years of data

% Rows are organized by:
% 1) prevalence (mean) of light infections (pre-SAC)
% 2) prevalence (SD) of light infections (pre-SAC)
% 3) prevalence (mean) of heavy infections (pre-SAC)
% 4) prevalence (SD) of heavy infections (pre-SAC)
% 5) prevalence (mean) of light infections (SAC)
% 6) prevalence (SD) of light infections (SAC)
% 7) prevalence (mean) of heavy infections (SAC)
% 8) prevalence (SD) of heavy infections (SAC)
% 9) prevalence (mean) of light infections (adult)
% 10) prevalence (SD) of light infections (adult)
% 11) prevalence (mean) of heavy infections (adult)
% 12) prevalence (SD) of heavy infections (adult)

% Excel files (folder: Stanford research -> Case Western -> Data)
% 1) PrvHighRist.xlsx
% 2) PrvLowRist.xlsx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% 1) total_DALY_sequelae_pre*- DALYs for each sequelae in pre-SAC
% 2) total_DALY_sequelae_c*- DALYs for each sequelae in SAC
% 3) total_DALY_sequelae_a*- DALYs for each sequelae in adults
% 4) ICER_values (1x2)- Total cost (USD) and total disability (DALY)
% 5) total_DALYs_age- Total DALYs for each age: 1) pre-SAC; 2) SAC; 3)
%    adults
% 6) total_cost- Total cost (USD)
% 7) total_DALY_iteration- DALY breakdown by time step (5 years, with
%    one month time step)
% *Not discounted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Case Western data inputs (based on tested strategies)
% MDA_strategy=[3, 3, 1, 1, 0, 0, 3, 3, 3, 3, 1, 1, 1, 1];
% snail_strategy=[0, 0, 0, 0, 1, 1, 1, 1 ,1 ,1 ,1 ,1 ,1 ,1];
% g_MDA=[0.75 0.75 0.75];
% g_snail=0.9;
% MDA_freq=[1, 2, 1, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 2];
% snail_freq=[0, 0, 0, 0, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2];

% 1.	CWA
% 2.	CWB
% 3.	SBA
% 4.	SBB
% 5.	SCA
% 6.	SCB
% 7.	CWA+SCA
% 8.	CWA+SCB
% 9.	CWB+SCA
% 10.	CWB+SCB
% 11.	SBA+SCA
% 12.	SBA+SCB
% 13.	SBB+SCA
% 14.	SBB+SCB



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

low_high=1; % (1)0 low prevalence- Mwaembe; (2) high prevalence- Milalina

lower_SD=0; % -1.96
upper_SD=0; % 1.96

% Non-compliance (0- missing at random; 0.1- 10% systematic noncompliance
% (base case scenario))
non_compliance=0.1;

% e) Discounting- 3% annual: (0) no; (1); yes
discounting_on=1; % keep on
disc_rate=0.03; %discounting 0.03

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% a) Costing inputs (2016 US$) and coverage

% Cost inputs (MDA)
cost_PRZ=0.21;
%cost_ALB=0.03;
base_cost=0.5;
community_cost_multiplier=3;

school_delivery_cost=base_cost;
community_delivery_cost=school_delivery_cost*community_cost_multiplier;
% community_delivery_cost=1.50;


% cost_children=cost_PRZ+cost_ALB+school_delivery_cost;
% cost_adults=cost_PRZ+cost_ALB+community_delivery_cost;
% cost_pre=cost_PRZ+cost_ALB+community_delivery_cost;
cost_children=cost_PRZ+school_delivery_cost;
cost_adults=cost_PRZ+community_delivery_cost;
cost_pre=cost_PRZ+community_delivery_cost;

% Cost inputs (molluscicides)
cost_snail_control_community=379.43;

% Coverage rate for mass treatment (g)
g_MDA_pre=g_MDA(1);
g_MDA_c=g_MDA(2);
g_MDA_a=g_MDA(3);

% Coverage rate for snail control
g_snail;

% b) Demographic and model parameters

% Proportion of population pre-SAC (X_pre), children (X_c), and adults (X_a)
% Case Western- Kenya stats
X_c=0.28; %0.153;
X_a=0.54 ; %0.763;
X_pre=0.18; % 0.084;

% Total population size for simulation
total_pop=5000;

% Compute size of each age population
total_pop_c=total_pop*X_c;
total_pop_a=total_pop*X_a;
total_pop_pre=total_pop*X_pre;

% c) Anemia modeling parameters

% Baseline Hb mean and SD for population in Cote d'Ivoire
mean_Hb_pop_c=112; SD_Hb_pop_c=15;
mean_Hb_pop_pre=mean_Hb_pop_c; SD_Hb_pop_pre=SD_Hb_pop_c;
mean_Hb_pop_a_men=134; SD_Hb_pop_a_men=19;
mean_Hb_pop_a_women=111; SD_Hb_pop_a_women=16;

% d) Caclulate time length of simulation
% see below

% f) Schistosomiasis disability and mortality rate
schisto_DALY= [0.014 0.02 0.05];
schisto_weight=schisto_DALY;
anemia_weight=[0.0041 0.0056 0.1615];

mortality_rate_schisto=1/20000; % Using 2010 GBD Lancet study on
% global mortality data (1 in 20,000 infections)

LE_CI=61; % Life expectancy in Kenya



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load transmission data
% Read Case Western data (excel)

if low_high==1 && non_compliance==0.1
    filename1 = 'Mwaembe 14.xls';
    
elseif low_high==2 && non_compliance==0.1
    filename1 = 'Milalani 14.xls';
    
elseif low_high==1 && non_compliance==0
    %filename1 = 'SWB projections of prevalences in LOW RISK Sh villages LDC population profile missing at random 091916 .xlsx';
    
elseif low_high==2 && non_compliance==0
    %filename1 = 'SWB projections of prevalences in HIGH RISK Sh villages LDC population profile missing at random 091916 .xlsx';
    
elseif low_high==2 && non_compliance==0.05
    %filename1 = 'SWB projections of prevalences in HIGH RISK Sh villages LDC population profile 5pc chronic non adherence 091916 .xlsx';
    
end


sheet_name = {'Strategy1', 'Strategy2', 'Strategy3', 'Strategy4', 'Strategy5', 'Strategy6', 'Strategy7', 'Strategy8', 'Strategy9', 'Strategy10', 'Strategy11', 'Strategy12', 'Strategy13', 'Strategy14'};

xlRange = 'A2:L122';

% Load transmission data into 3D array (12x121x14)
% Note: 3rd dimension (z) is for each of 14 strategies

for i=1:length(sheet_name)
    Data_High_temp(:,:,i)= xlsread(filename1, char(sheet_name(i)),xlRange);
end

Data_High_1 = zeros(size(Data_High_temp)+[0,0,1]);

Data_High_1(:,:,1) = repmat(Data_High_temp(1,:,1),121,1);
Data_High_1(:,:,2:15) = Data_High_temp;


% d) Caclulate time length of simulation
epi_vector_size=length(Data_High_1);
time_model=ceil(epi_vector_size(1));

% Re-format data for CEA code.
% Rows: Monthly measurements
% Columns: 1) prevalence (0-1); 2) number uninfected (#); 3) number light
% infections (#); 4) number heavy infections (#)

% Mean of transmission results
for i=1:size(Data_High_1, 3)
    
    % Excel_decision
    % 1- which sheet
    % 2- base case (0), lower (-1), upper(1)
    schisto_epi_vector_pre_full(:,3:4, i) = (Data_High_1(:,[1,7],i));
    schisto_epi_vector_children_full(:,3:4, i) = (Data_High_1(:,[3,9],i));
    schisto_epi_vector_adults_full(:,3:4, i) = (Data_High_1(:,[5,11],i));
    schisto_epi_SD_pre= (Data_High_1(:,[2,8],i));
    schisto_epi_SD_c= (Data_High_1(:,[4,10],i));
    schisto_epi_SD_a= (Data_High_1(:,[6,12],i));
    
    if excel_decision(2)==0
        
    elseif excel_decision(2)==-1
        schisto_epi_vector_pre_full(:,3:4, i)=schisto_epi_vector_pre_full(:,3:4, i)- lower_SD*schisto_epi_SD_pre;
        schisto_epi_vector_children_full(:,3:4, i)=schisto_epi_vector_children_full(:,3:4, i)- lower_SD*schisto_epi_SD_c;
        schisto_epi_vector_adults_full(:,3:4, i)=schisto_epi_vector_adults_full(:,3:4, i)- lower_SD*schisto_epi_SD_a;
        
    elseif excel_decision(2)==1
        schisto_epi_vector_pre_full(:,3:4, i)=schisto_epi_vector_pre_full(:,3:4, i)+ upper_SD*schisto_epi_SD_pre;
        schisto_epi_vector_children_full(:,3:4, i)=schisto_epi_vector_children_full(:,3:4, i)+ upper_SD*schisto_epi_SD_c;
        schisto_epi_vector_adults_full(:,3:4, i)=schisto_epi_vector_adults_full(:,3:4, i)+ upper_SD*schisto_epi_SD_a;
        
    end
    
    
    schisto_epi_vector_pre_full(:,1, i)= schisto_epi_vector_pre_full(:,3, i) + schisto_epi_vector_pre_full(:,4, i);
    schisto_epi_vector_pre_full(:,2, i)= round((1-schisto_epi_vector_pre_full(:,1, i))*total_pop_pre);
    schisto_epi_vector_pre_full(:,3, i)= round(schisto_epi_vector_pre_full(:,3, i)*total_pop_pre);
    schisto_epi_vector_pre_full(:,4, i)= round(schisto_epi_vector_pre_full(:,4, i)*total_pop_pre);
    
    schisto_epi_vector_children_full(:,1, i)= schisto_epi_vector_children_full(:,3, i) + schisto_epi_vector_children_full(:,4, i);
    schisto_epi_vector_children_full(:,2, i)= round((1-schisto_epi_vector_children_full(:,1, i))*total_pop_c);
    schisto_epi_vector_children_full(:,3, i)= round(schisto_epi_vector_children_full(:,3, i)*total_pop_c);
    schisto_epi_vector_children_full(:,4, i)= round(schisto_epi_vector_children_full(:,4, i)*total_pop_c);
    
    schisto_epi_vector_adults_full(:,1, i)= schisto_epi_vector_adults_full(:,3, i) + schisto_epi_vector_adults_full(:,4, i);
    schisto_epi_vector_adults_full(:,2, i)= round((1-schisto_epi_vector_adults_full(:,1, i))*total_pop_a);
    schisto_epi_vector_adults_full(:,3, i)= round(schisto_epi_vector_adults_full(:,3, i)*total_pop_a);
    schisto_epi_vector_adults_full(:,4, i)= round(schisto_epi_vector_adults_full(:,4, i)*total_pop_a);
    
    % Sensitivity analysis: Can take mean+/-SD of transmission results to
    % generate uncertainty interval
end

schisto_epi_vector_pre = schisto_epi_vector_pre_full(:,:,excel_decision(1));
schisto_epi_vector_children = schisto_epi_vector_children_full(:,:,excel_decision(1));
schisto_epi_vector_adults = schisto_epi_vector_adults_full(:,:,excel_decision(1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute mortality

% Compute average number of years of life lost in case of death
% resulting from schistosomiasis. An average discounted DALY value is
% assigned based upon whether death occurs in pre-SAC, SAC, or adult
% population.
DALY_schisto_death_a_years=round(LE_CI-((15+LE_CI)/2));
DALY_schisto_death_c_years=round(LE_CI-((5+14)/2));
DALY_schisto_death_pre_years=round(LE_CI-((1+4)/2));

% Years of life (DALYs lost) is discounted at 3% annually as per
% convention.
for t=1:LE_CI
    if discounting_on==1; discount=1/((1+(disc_rate))^(t-1)); end
    if discounting_on==0; discount=1; end
    if t<=DALY_schisto_death_a_years
        DALY_schisto_death_a(t)=1*discount;
    end
    
    if t<=DALY_schisto_death_c_years
        DALY_schisto_death_c(t)=1*discount;
    end
    
    if t<=DALY_schisto_death_pre_years
        DALY_schisto_death_pre(t)=1*discount;
    end
end

DALY_schisto_death_a=sum(DALY_schisto_death_a);
DALY_schisto_death_c=sum(DALY_schisto_death_c);
DALY_schisto_death_pre=sum(DALY_schisto_death_pre);

drawnow; pause(0.1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Disability modeling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Disability #1: Anemia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main Steps:
% 1) Input parameters
% 2) Mixture model methodology
% 3) Quantify baseline anemia
% 4) Compute anemia at each time step

% 1) Input parameters
% Anemia thresold (Hb; g/L)
% Anemia thresholds for pre-SAC
pre_anemia_low=110;
pre_anemia_mod=100;
pre_anemia_severe=70;

% Anemia thresholds for SAC
SAC_anemia_low=116.5;
SAC_anemia_mod=110;
SAC_anemia_severe=80;

% Anemia thresholds for men
men_anemia_low=130;
men_anemia_mod=110;
men_anemia_severe=80;

% Anemia thresholds for women
women_anemia_low=120;
women_anemia_mod=110;
women_anemia_severe=80;

PZQ_eff=2.8;
%PZQ_eff=0;

PZQ_treatment_effect_a=PZQ_eff;
PZQ_treatment_effect_c=PZQ_eff;
PZQ_treatment_effect_pre=PZQ_eff;

% 2) Mixture model methodology
% Goal is to to discern Hb of uninfected individuals for all
% sub-populations at baseline

% 1) SAC
% Calculate proportion of children who are: 1) Schisto infected
baseline_frac_c(1,1)= schisto_epi_vector_children(1,1);
% Perform mixture model calculations to discern Hb of uninfected children
% for future calculations
% The following set of two equations are used in this calculation:
% 1) Hb_uninfected-Hb_schisto=2.80
% 3) Hb_uninfected*X_uninfected + Hb_schisto*X_schisto= Obs_pop_Hb
Hb_matrix_baseline_c=[1 -1;  (1-baseline_frac_c(1,1)) baseline_frac_c(1,1)];
Hb_Bmatrix_baseline_c=[2.80 mean_Hb_pop_c]';
Hb_matrix_results_c=(Hb_matrix_baseline_c\Hb_Bmatrix_baseline_c)';
Hb_uninfected_c=Hb_matrix_results_c(1);

% 2) Pre-SAC
% Calculate proportion of preSAC who are: 1) Schisto infected
baseline_frac_pre(1,1)= schisto_epi_vector_pre(1,1);
% Perform mixture model calculations to discern Hb of uninfected preSAC
% for future calculations
% The following set of two equations are used in this calculation:
% 1) Hb_uninfected-Hb_schisto=2.80
% 3) Hb_uninfected*X_uninfected + Hb_schisto*X_schisto= Obs_pop_Hb
Hb_matrix_baseline_pre=[1 -1; (1-baseline_frac_pre(1,1)) baseline_frac_pre(1,1)];
Hb_Bmatrix_baseline_pre=[2.80 mean_Hb_pop_pre]';
Hb_matrix_results_pre=(Hb_matrix_baseline_pre\Hb_Bmatrix_baseline_pre)';
Hb_uninfected_pre=Hb_matrix_results_pre(1);


% 3) Adult- males
% Calculate proportion of men who are: 1) Schisto infected
baseline_frac_a_men(1,1)= schisto_epi_vector_adults(1,1);
% Perform mixture model calculations to discern Hb of uninfected men
% for future calculations
% The following set of two equations are used in this calculation:
% 1) Hb_uninfected-Hb_schisto=2.80
% 3) Hb_uninfected*X_uninfected + Hb_schisto*X_schisto= Obs_pop_Hb
Hb_matrix_baseline_a_men=[1 -1; (1-baseline_frac_a_men(1,1)) baseline_frac_a_men(1,1)];
Hb_Bmatrix_baseline_a_men=[2.80 mean_Hb_pop_a_men]';
Hb_matrix_results_a_men=(Hb_matrix_baseline_a_men\Hb_Bmatrix_baseline_a_men)';
Hb_uninfected_a_men=Hb_matrix_results_a_men(1);

% 4) Adult- females
% Calculate proportion of men who are: 1) Schisto infected
baseline_frac_a_women(1,1)= schisto_epi_vector_adults(1,1);
% Perform mixture model calculations to discern Hb of uninfected women
% for future calculations
% The following set of two equations are used in this calculation:
% 1) Hb_uninfected-Hb_schisto=2.80
% 3) Hb_uninfected*X_uninfected + Hb_schisto*X_schisto= Obs_pop_Hb
Hb_matrix_baseline_a_women=[1 -1; (1-baseline_frac_a_women(1,1)) baseline_frac_a_women(1,1)];
Hb_Bmatrix_baseline_a_women=[2.80 mean_Hb_pop_a_women]';
Hb_matrix_results_a_women=(Hb_matrix_baseline_a_women\Hb_Bmatrix_baseline_a_women)';
Hb_uninfected_a_women=Hb_matrix_results_a_women(1);


% Update mean_Hb_pop to reflect mean Hb of uninfected individuals
mean_Hb_pop_c=Hb_uninfected_c;
mean_Hb_pop_a_men=Hb_uninfected_a_men;
mean_Hb_pop_a_women=Hb_uninfected_a_women;
mean_Hb_pop_pre=Hb_uninfected_pre;

% 3) Quantify baseline anemia in worm free population
% Quantify baseline anemia (i.e. in a population without schisto/helminths, what is the baseline anemia)
% In these computations, a schisto-free mean Hb is used for each population
%

% Matrix organization
% 1- severe anemia
% 2- moderate anemia
% 3- mild anemia
% 4- no anemia

% Children
% Heavy anemia
anemia_tracker_c_baseline(1,1)=round(normcdf(SAC_anemia_severe,mean_Hb_pop_c, SD_Hb_pop_c)*total_pop_c);
% Moderate anemia
anemia_tracker_c_baseline(1,2)=round(normcdf(SAC_anemia_mod,mean_Hb_pop_c, SD_Hb_pop_c)*total_pop_c)-anemia_tracker_c_baseline(1,1);
% Heavy anemia
anemia_tracker_c_baseline(1,3)=round(normcdf(SAC_anemia_low,mean_Hb_pop_c, SD_Hb_pop_c)*total_pop_c)-anemia_tracker_c_baseline(1,1)-anemia_tracker_c_baseline(1,2);
% No anemia
anemia_tracker_c_baseline(1,4)=round(normcdf(inf,mean_Hb_pop_c, SD_Hb_pop_c)*total_pop_c)-anemia_tracker_c_baseline(1,1)-anemia_tracker_c_baseline(1,2)-anemia_tracker_c_baseline(1,3);

% Pre-SAC
% Heavy anemia
anemia_tracker_pre_baseline(1,1)=round(normcdf(pre_anemia_severe,mean_Hb_pop_pre, SD_Hb_pop_c)*total_pop_pre);
% Moderate anemia
anemia_tracker_pre_baseline(1,2)=round(normcdf(pre_anemia_mod,mean_Hb_pop_pre, SD_Hb_pop_c)*total_pop_pre)-anemia_tracker_pre_baseline(1,1);
% Heavy anemia
anemia_tracker_pre_baseline(1,3)=round(normcdf(pre_anemia_low,mean_Hb_pop_pre, SD_Hb_pop_c)*total_pop_pre)-anemia_tracker_pre_baseline(1,1)-anemia_tracker_pre_baseline(1,2);
% No anemia
anemia_tracker_pre_baseline(1,4)=round(normcdf(inf,mean_Hb_pop_pre, SD_Hb_pop_c)*total_pop_pre)-anemia_tracker_pre_baseline(1,1)-anemia_tracker_pre_baseline(1,2)-anemia_tracker_pre_baseline(1,3);

% Adults-male
% Heavy anemia
anemia_tracker_a_men_baseline(1,1)=round(normcdf(men_anemia_severe,mean_Hb_pop_a_men, SD_Hb_pop_a_men)*total_pop_a*0.5);
% Moderate anemia
anemia_tracker_a_men_baseline(1,2)=round(normcdf(men_anemia_mod,mean_Hb_pop_a_men, SD_Hb_pop_a_men)*total_pop_a*0.5)-anemia_tracker_a_men_baseline(1,1);
% Heavy anemia
anemia_tracker_a_men_baseline(1,3)=round(normcdf(men_anemia_low,mean_Hb_pop_a_men, SD_Hb_pop_a_men)*total_pop_a*0.5)-anemia_tracker_a_men_baseline(1,1)-anemia_tracker_a_men_baseline(1,2);
% No anemia
anemia_tracker_a_men_baseline(1,4)=round(normcdf(inf,mean_Hb_pop_a_men, SD_Hb_pop_a_men)*total_pop_a*0.5)-anemia_tracker_a_men_baseline(1,1)-anemia_tracker_a_men_baseline(1,2)-anemia_tracker_a_men_baseline(1,3);

% Adults-female
% Heavy anemia
anemia_tracker_a_women_baseline(1,1)=round(normcdf(women_anemia_severe,mean_Hb_pop_a_women, SD_Hb_pop_a_women)*total_pop_a*0.5);
% Moderate anemia
anemia_tracker_a_women_baseline(1,2)=round(normcdf(women_anemia_mod,mean_Hb_pop_a_women, SD_Hb_pop_a_women)*total_pop_a*0.5)-anemia_tracker_a_women_baseline(1,1);
% Heavy anemia
anemia_tracker_a_women_baseline(1,3)=round(normcdf(women_anemia_low,mean_Hb_pop_a_women, SD_Hb_pop_a_women)*total_pop_a*0.5)-anemia_tracker_a_women_baseline(1,1)-anemia_tracker_a_women_baseline(1,2);
% No anemia
anemia_tracker_a_women_baseline(1,4)=round(normcdf(inf,mean_Hb_pop_a_women, SD_Hb_pop_a_women)*total_pop_a*0.5)-anemia_tracker_a_women_baseline(1,1)-anemia_tracker_a_women_baseline(1,2)-anemia_tracker_a_women_baseline(1,3);

% Adults- combined
anemia_tracker_a_combined_baseline(1,:)=anemia_tracker_a_men_baseline(1,:)+anemia_tracker_a_women_baseline(1,:);

% 4) Compute anemia at each time step
for i= 1:time_model
    
    % Calculate anemia present at each time step of the 5 year simulation,
    % including light, moderate, and heavy anemia following
    % age/gender-specific hemoglobin (Hb) thresholds as previously stated.
    % First, the proportion of each age groups that is: (1) uninfected; (2)
    % infected w/ schisto; (3) the mean Hb is calculated based upon
    % previously calculated "uninfected" Hb, a  schisto Hb downshift of
    % 2.80 g/L following the GBD anemia study (Blood, 2014)
    
    % SAC
    anemia_matrix_c(i,1)= 1-schisto_epi_vector_children(i,1);
    anemia_matrix_c(i,2)= schisto_epi_vector_children(i,1);
    anemia_matrix_c(i,4)= mean_Hb_pop_c*anemia_matrix_c(i,1) + (mean_Hb_pop_c-PZQ_treatment_effect_c)*anemia_matrix_c(i,2);
    
    % Pre-SAC
    anemia_matrix_pre(i,1)= 1-schisto_epi_vector_pre(i,1);
    anemia_matrix_pre(i,2)= schisto_epi_vector_pre(i,1);
    anemia_matrix_pre(i,4)= mean_Hb_pop_pre*anemia_matrix_pre(i,1) + (mean_Hb_pop_pre-PZQ_treatment_effect_pre)*anemia_matrix_pre(i,2);
    
    % Adults- Male
    anemia_matrix_a_men(i,1)= 1-schisto_epi_vector_adults(i,1);
    anemia_matrix_a_men(i,2)= schisto_epi_vector_adults(i,1);
    anemia_matrix_a_men(i,4)= mean_Hb_pop_a_men*anemia_matrix_a_men(i,1) + (mean_Hb_pop_a_men-PZQ_treatment_effect_a)*anemia_matrix_a_men(i,2);
    
    % Adults- Female
    anemia_matrix_a_women(i,1)= 1-schisto_epi_vector_adults(i,1);
    anemia_matrix_a_women(i,2)= schisto_epi_vector_adults(i,1);
    anemia_matrix_a_women(i,4)= mean_Hb_pop_a_women*anemia_matrix_a_women(i,1) + (mean_Hb_pop_a_women-PZQ_treatment_effect_a)*anemia_matrix_a_women(i,2);
    
    % Calculate anemia at each iteration
    % 1- heavy anemia
    % 2- moderate anemia
    % 3- light anemia
    % 4- no anemia
    
    % a) SAC
    % Heavy anemia
    anemia_tracker_c(i,1)=round(normcdf(SAC_anemia_severe,anemia_matrix_c(i,4), SD_Hb_pop_c)*total_pop_c);
    % Moderate anemia
    anemia_tracker_c(i,2)=round(normcdf(SAC_anemia_mod,anemia_matrix_c(i,4), SD_Hb_pop_c)*total_pop_c)-anemia_tracker_c(i,1);
    % Light anemia
    anemia_tracker_c(i,3)=round(normcdf(SAC_anemia_low,anemia_matrix_c(i,4), SD_Hb_pop_c)*total_pop_c)-anemia_tracker_c(i,1)-anemia_tracker_c(i,2);
    % No anemia
    anemia_tracker_c(i,4)=round(normcdf(inf,anemia_matrix_c(i,4), SD_Hb_pop_c)*total_pop_c)-anemia_tracker_c(i,1)-anemia_tracker_c(i,2)-anemia_tracker_c(i,3);
    
    % b) Pre-SAC
    % Heavy anemia
    anemia_tracker_pre(i,1)=round(normcdf(pre_anemia_severe,anemia_matrix_pre(i,4), SD_Hb_pop_pre)*total_pop_pre);
    % Moderate anemia
    anemia_tracker_pre(i,2)=round(normcdf(pre_anemia_mod,anemia_matrix_pre(i,4), SD_Hb_pop_pre)*total_pop_pre)-anemia_tracker_pre(i,1);
    % Light anemia
    anemia_tracker_pre(i,3)=round(normcdf(pre_anemia_low,anemia_matrix_pre(i,4), SD_Hb_pop_pre)*total_pop_pre)-anemia_tracker_pre(i,1)-anemia_tracker_pre(i,2);
    % No anemia
    anemia_tracker_pre(i,4)=round(normcdf(inf,anemia_matrix_pre(i,4), SD_Hb_pop_pre)*total_pop_pre)-anemia_tracker_pre(i,1)-anemia_tracker_pre(i,2)-anemia_tracker_pre(i,3);
    
    % c) Adults-male
    % Heavy anemia
    anemia_tracker_a_men(i,1)=round(normcdf(men_anemia_severe,anemia_matrix_a_men(i,4), SD_Hb_pop_a_men)*total_pop_a*0.5);
    % Moderate anemia
    anemia_tracker_a_men(i,2)=round(normcdf(men_anemia_mod,anemia_matrix_a_men(i,4), SD_Hb_pop_a_men)*total_pop_a*0.5)-anemia_tracker_a_men(i,1);
    % Light anemia
    anemia_tracker_a_men(i,3)=round(normcdf(men_anemia_low,anemia_matrix_a_men(i,4), SD_Hb_pop_a_men)*total_pop_a*0.5)-anemia_tracker_a_men(i,1)-anemia_tracker_a_men(i,2);
    % No anemia
    anemia_tracker_a_men(i,4)=round(normcdf(inf,anemia_matrix_a_men(i,4), SD_Hb_pop_a_men)*total_pop_a*0.5)-anemia_tracker_a_men(i,1)-anemia_tracker_a_men(i,2)-anemia_tracker_a_men(i,3);
    
    % d) Adults-female
    % Heavy anemia
    anemia_tracker_a_women(i,1)=round(normcdf(women_anemia_severe,anemia_matrix_a_women(i,4), SD_Hb_pop_a_women)*total_pop_a*0.5);
    % Moderate anemia
    anemia_tracker_a_women(i,2)=round(normcdf(women_anemia_mod,anemia_matrix_a_women(i,4), SD_Hb_pop_a_women)*total_pop_a*0.5)-anemia_tracker_a_women(i,1);
    % Light anemia
    anemia_tracker_a_women(i,3)=round(normcdf(women_anemia_low,anemia_matrix_a_women(i,4), SD_Hb_pop_a_women)*total_pop_a*0.5)-anemia_tracker_a_women(i,1)-anemia_tracker_a_women(i,2);
    % No anemia
    anemia_tracker_a_women(i,4)=round(normcdf(inf,anemia_matrix_a_women(i,4), SD_Hb_pop_a_women)*total_pop_a*0.5)-anemia_tracker_a_women(i,1)-anemia_tracker_a_women(i,2)-anemia_tracker_a_women(i,3);
    
    % Adults- combined
    anemia_tracker_a_combined(i,:)=anemia_tracker_a_men(i,:)+anemia_tracker_a_women(i,:);
    
end

for i=2:time_model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Disability #1b: Anemia counting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Counter anemia is a matrix that stores averted anemia: 1- severe anemia
    % (0.1615 DALYs); 2- moderate anemia (0.0056 DALYs), 3- light anemia
    % (0.0041 DALYs)
    
    % Severe anemia counter
    counter_anemia_adults(i,1)= anemia_tracker_a_combined(i,1)-anemia_tracker_a_combined_baseline(1,1);
    counter_anemia_children(i,1)= anemia_tracker_c(i,1)-anemia_tracker_c_baseline(1,1);
    counter_anemia_pre(i,1)= anemia_tracker_pre(i,1)-anemia_tracker_pre_baseline(1,1);
    
    % Moderate anemia counter
    counter_anemia_adults(i,2)= anemia_tracker_a_combined(i,2)-anemia_tracker_a_combined_baseline(1,2);
    counter_anemia_children(i,2)= anemia_tracker_c(i,2)-anemia_tracker_c_baseline(1,2);
    counter_anemia_pre(i,2)= anemia_tracker_pre(i,2)-anemia_tracker_pre_baseline(1,2);
    
    % Light anemia counter
    counter_anemia_adults(i,3)= anemia_tracker_a_combined(i,3)-anemia_tracker_a_combined_baseline(1,3);
    counter_anemia_children(i,3)= anemia_tracker_c(i,3)-anemia_tracker_c_baseline(1,3);
    counter_anemia_pre(i,3)= anemia_tracker_pre(i,3)-anemia_tracker_pre_baseline(1,3);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Disability #2: Schistosomiasis (S. haematobium) symptomatic infection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DALY weight applied differentially based upon infection intensity
    % light- 1.4% and heavy- 5%
    
    schisto_sequela_c(i,1) =  (schisto_epi_vector_children(i,3));
    schisto_sequela_a(i,1) =  (schisto_epi_vector_adults(i,3));
    schisto_sequela_pre(i,1)= (schisto_epi_vector_pre(i,3));
    
    schisto_sequela_c(i,2) =  (schisto_epi_vector_children(i,4));
    schisto_sequela_a(i,2) =  (schisto_epi_vector_adults(i,4));
    schisto_sequela_pre(i,2)= (schisto_epi_vector_pre(i,4));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Distribution of disease within population
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Creates matrix with DALY weight assigned in a row corresponding to the
    % probabilistic combination of sequelae. Deterministic
    DALY_matrix_count_c=zeros(1, 10);
    DALY_matrix_count_a=zeros(1, 10);
    DALY_matrix_count_pre=zeros(1, 10);
    
    % 1- schistosomiasis
    % 2- Anemia (three transitions possible)
    % 3- total DALYs
    
    % Calculates probability of each sequelae.
    
    % Note that probabilities are computed over denominator of all
    % people
    schisto_prob1_c= schisto_sequela_c(i,1)/total_pop_c;
    schisto_prob1_a= schisto_sequela_a(i,1)/total_pop_a;
    schisto_prob1_pre= schisto_sequela_pre(i,1)/total_pop_pre;
    
    schisto_prob2_c= schisto_sequela_c(i,2)/total_pop_c;
    schisto_prob2_a= schisto_sequela_a(i,2)/total_pop_a;
    schisto_prob2_pre= schisto_sequela_pre(i,2)/total_pop_pre;
    
    schisto_deathprob1_c= schisto_epi_vector_children(i,1);
    schisto_deathprob1_a= schisto_epi_vector_adults(i,1);
    schisto_deathprob1_pre= schisto_epi_vector_pre(i,1);
    
    % Note that probabilities are computed over denominator of infected
    % people
    prob_anemia_low_a=counter_anemia_adults(i,3)/(total_pop_a*schisto_epi_vector_adults(i,1));
    prob_anemia_low_c=counter_anemia_children(i,3)/(total_pop_c*schisto_epi_vector_children(i,1));
    prob_anemia_low_pre=counter_anemia_pre(i,3)/(total_pop_pre*schisto_epi_vector_pre(i,1));
    
    prob_anemia_mod_a=counter_anemia_adults(i,2)/(total_pop_a*schisto_epi_vector_adults(i,1));
    prob_anemia_mod_c=counter_anemia_children(i,2)/(total_pop_c*schisto_epi_vector_children(i,1));
    prob_anemia_mod_pre=counter_anemia_pre(i,2)/(total_pop_pre*schisto_epi_vector_pre(i,1));
    
    prob_anemia_severe_a=counter_anemia_adults(i,1)/(total_pop_a*schisto_epi_vector_adults(i,1));
    prob_anemia_severe_c=counter_anemia_children(i,1)/(total_pop_c*schisto_epi_vector_children(i,1));
    prob_anemia_severe_pre=counter_anemia_pre(i,1)/(total_pop_pre*schisto_epi_vector_pre(i,1));
    
    
    
    % Define disability (multiplicative)
    disability_schisto_no_anemia_no= 0;
    disability_schisto_light_anemia_no= 1 - ( (1-schisto_weight(1))*(1-0));
    disability_schisto_light_anemia_low= 1 - ( (1-schisto_weight(1))*(1-anemia_weight(1)));
    disability_schisto_light_anemia_mod= 1 - ( (1-schisto_weight(1))*(1-anemia_weight(2)));
    disability_schisto_light_anemia_sev= 1 - ( (1-schisto_weight(1))*(1-anemia_weight(3)));
    disability_schisto_heavy_anemia_no= 1 - ( (1-schisto_weight(3))*(1-0));
    disability_schisto_heavy_anemia_low= 1 - ( (1-schisto_weight(3))*(1-anemia_weight(1)));
    disability_schisto_heavy_anemia_mod= 1 - ( (1-schisto_weight(3))*(1-anemia_weight(2)));
    disability_schisto_heavy_anemia_sev= 1 - ( (1-schisto_weight(3))*(1-anemia_weight(3)));
    
    % Assign probablity of each sequelae
    % 1) Assign each sequelae for ADULTS
    % Prob (schisto) (anemia)
    prob_schisto_no_anemia_no_a= (1-schisto_prob1_a-schisto_prob2_a);
    prob_schisto_light_anemia_no_a= schisto_prob1_a * (1-prob_anemia_low_a-prob_anemia_mod_a-prob_anemia_severe_a);
    prob_schisto_light_anemia_low_a= schisto_prob1_a * (prob_anemia_low_a);
    prob_schisto_light_anemia_mod_a= schisto_prob1_a * (prob_anemia_mod_a);
    prob_schisto_light_anemia_sev_a= schisto_prob1_a * (prob_anemia_severe_a);
    prob_schisto_heavy_anemia_no_a= schisto_prob2_a * (1-prob_anemia_low_a-prob_anemia_mod_a-prob_anemia_severe_a);
    prob_schisto_heavy_anemia_low_a= schisto_prob2_a * (prob_anemia_low_a);
    prob_schisto_heavy_anemia_mod_a= schisto_prob2_a * (prob_anemia_mod_a);
    prob_schisto_heavy_anemia_sev_a= schisto_prob2_a * (prob_anemia_severe_a);
    
    prob_total_a= (prob_schisto_no_anemia_no_a + prob_schisto_light_anemia_no_a + prob_schisto_light_anemia_low_a + prob_schisto_light_anemia_mod_a + prob_schisto_light_anemia_sev_a + prob_schisto_heavy_anemia_no_a + prob_schisto_heavy_anemia_low_a +  prob_schisto_heavy_anemia_mod_a + prob_schisto_heavy_anemia_sev_a);
    
    % Disability from each combination of sequelae (including no
    % disability). Multiplicative accounted for multiple sequelae.
    DALY_matrix_count_a(1)= prob_schisto_no_anemia_no_a * total_pop_a * disability_schisto_no_anemia_no;
    DALY_matrix_count_a(2)= prob_schisto_light_anemia_no_a * total_pop_a * disability_schisto_light_anemia_no;
    DALY_matrix_count_a(3)= prob_schisto_light_anemia_low_a * total_pop_a * disability_schisto_light_anemia_low;
    DALY_matrix_count_a(4)= prob_schisto_light_anemia_mod_a * total_pop_a * disability_schisto_light_anemia_mod;
    DALY_matrix_count_a(5)= prob_schisto_light_anemia_sev_a * total_pop_a * disability_schisto_light_anemia_sev;
    DALY_matrix_count_a(6)= prob_schisto_heavy_anemia_no_a * total_pop_a * disability_schisto_heavy_anemia_no;
    DALY_matrix_count_a(7)= prob_schisto_heavy_anemia_low_a * total_pop_a * disability_schisto_heavy_anemia_low;
    DALY_matrix_count_a(8)= prob_schisto_heavy_anemia_mod_a * total_pop_a * disability_schisto_heavy_anemia_mod;
    DALY_matrix_count_a(9)= prob_schisto_heavy_anemia_sev_a * total_pop_a * disability_schisto_heavy_anemia_sev;
    
    % Compute mortality
    DALY_matrix_count_a(10)=(DALY_schisto_death_a*schisto_deathprob1_a*total_pop_a*mortality_rate_schisto);
    
    
    
    % 2) Assign each sequelae for CHILDREN
    % Prob (schisto) (anemia)
    prob_schisto_no_anemia_no_c= (1-schisto_prob1_c-schisto_prob2_c);
    prob_schisto_light_anemia_no_c= schisto_prob1_c * (1-prob_anemia_low_c-prob_anemia_mod_c-prob_anemia_severe_c);
    prob_schisto_light_anemia_low_c= schisto_prob1_c * (prob_anemia_low_c);
    prob_schisto_light_anemia_mod_c= schisto_prob1_c * (prob_anemia_mod_c);
    prob_schisto_light_anemia_sev_c= schisto_prob1_c * (prob_anemia_severe_c);
    prob_schisto_heavy_anemia_no_c= schisto_prob2_c * (1-prob_anemia_low_c-prob_anemia_mod_c-prob_anemia_severe_c);
    prob_schisto_heavy_anemia_low_c= schisto_prob2_c * (prob_anemia_low_c);
    prob_schisto_heavy_anemia_mod_c= schisto_prob2_c * (prob_anemia_mod_c);
    prob_schisto_heavy_anemia_sev_c= schisto_prob2_c * (prob_anemia_severe_c);
    
    prob_total_c= (prob_schisto_no_anemia_no_c + prob_schisto_light_anemia_no_c + prob_schisto_light_anemia_low_c + prob_schisto_light_anemia_mod_c + prob_schisto_light_anemia_sev_c + prob_schisto_heavy_anemia_no_c + prob_schisto_heavy_anemia_low_c +  prob_schisto_heavy_anemia_mod_c + prob_schisto_heavy_anemia_sev_c);
    
    % Disability from each combination of sequelae (including no
    % disability). Multiplicative accounted for multiple sequelae.
    DALY_matrix_count_c(1)= prob_schisto_no_anemia_no_c * total_pop_c * disability_schisto_no_anemia_no;
    DALY_matrix_count_c(2)= prob_schisto_light_anemia_no_c * total_pop_c * disability_schisto_light_anemia_no;
    DALY_matrix_count_c(3)= prob_schisto_light_anemia_low_c * total_pop_c * disability_schisto_light_anemia_low;
    DALY_matrix_count_c(4)= prob_schisto_light_anemia_mod_c * total_pop_c * disability_schisto_light_anemia_mod;
    DALY_matrix_count_c(5)= prob_schisto_light_anemia_sev_c * total_pop_c * disability_schisto_light_anemia_sev;
    DALY_matrix_count_c(6)= prob_schisto_heavy_anemia_no_c * total_pop_c * disability_schisto_heavy_anemia_no;
    DALY_matrix_count_c(7)= prob_schisto_heavy_anemia_low_c * total_pop_c * disability_schisto_heavy_anemia_low;
    DALY_matrix_count_c(8)= prob_schisto_heavy_anemia_mod_c * total_pop_c * disability_schisto_heavy_anemia_mod;
    DALY_matrix_count_c(9)= prob_schisto_heavy_anemia_sev_c * total_pop_c * disability_schisto_heavy_anemia_sev;
    
    % Compute mortality
    DALY_matrix_count_c(10)=(DALY_schisto_death_c*schisto_deathprob1_c*total_pop_c*mortality_rate_schisto);
    
    
    
    % 3) Assign probablity of each sequelae for Pre-SAC
    % Prob (schisto) (anemia)
    prob_schisto_no_anemia_no_pre= (1-schisto_prob1_pre-schisto_prob2_pre);
    prob_schisto_light_anemia_no_pre= schisto_prob1_pre * (1-prob_anemia_low_pre-prob_anemia_mod_pre-prob_anemia_severe_pre);
    prob_schisto_light_anemia_low_pre= schisto_prob1_pre * (prob_anemia_low_pre);
    prob_schisto_light_anemia_mod_pre= schisto_prob1_pre * (prob_anemia_mod_pre);
    prob_schisto_light_anemia_sev_pre= schisto_prob1_pre * (prob_anemia_severe_pre);
    prob_schisto_heavy_anemia_no_pre= schisto_prob2_pre * (1-prob_anemia_low_pre-prob_anemia_mod_pre-prob_anemia_severe_pre);
    prob_schisto_heavy_anemia_low_pre= schisto_prob2_pre * (prob_anemia_low_pre);
    prob_schisto_heavy_anemia_mod_pre= schisto_prob2_pre * (prob_anemia_mod_pre);
    prob_schisto_heavy_anemia_sev_pre= schisto_prob2_pre * (prob_anemia_severe_pre);
    
    prob_total_pre= (prob_schisto_no_anemia_no_pre + prob_schisto_light_anemia_no_pre + prob_schisto_light_anemia_low_pre + prob_schisto_light_anemia_mod_pre + prob_schisto_light_anemia_sev_pre + prob_schisto_heavy_anemia_no_pre + prob_schisto_heavy_anemia_low_pre +  prob_schisto_heavy_anemia_mod_pre + prob_schisto_heavy_anemia_sev_pre);
    
    % Disability from each combination of sequelae (including no
    % disability). Multiplicative accounted for multiple sequelae.
    DALY_matrix_count_pre(1)= prob_schisto_no_anemia_no_pre * total_pop_pre * disability_schisto_no_anemia_no;
    DALY_matrix_count_pre(2)= prob_schisto_light_anemia_no_pre * total_pop_pre * disability_schisto_light_anemia_no;
    DALY_matrix_count_pre(3)= prob_schisto_light_anemia_low_pre * total_pop_pre * disability_schisto_light_anemia_low;
    DALY_matrix_count_pre(4)= prob_schisto_light_anemia_mod_pre * total_pop_pre * disability_schisto_light_anemia_mod;
    DALY_matrix_count_pre(5)= prob_schisto_light_anemia_sev_pre * total_pop_pre * disability_schisto_light_anemia_sev;
    DALY_matrix_count_pre(6)= prob_schisto_heavy_anemia_no_pre * total_pop_pre * disability_schisto_heavy_anemia_no;
    DALY_matrix_count_pre(7)= prob_schisto_heavy_anemia_low_pre * total_pop_pre * disability_schisto_heavy_anemia_low;
    DALY_matrix_preount_pre(8)= prob_schisto_heavy_anemia_mod_pre * total_pop_pre * disability_schisto_heavy_anemia_mod;
    DALY_matrix_preount_pre(9)= prob_schisto_heavy_anemia_sev_pre * total_pop_pre * disability_schisto_heavy_anemia_sev;
    
    % Compute mortality
    DALY_matrix_count_pre(10)=(DALY_schisto_death_pre*schisto_deathprob1_pre*total_pop_pre*mortality_rate_schisto);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cost-effectiveness analysis: Disability
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Sum total DALYs per time step (month)
    total_DALY_iteration(i,:)=  [sum(DALY_matrix_count_pre)*(1/12),  sum(DALY_matrix_count_c)*(1/12), sum(DALY_matrix_count_a)*(1/12)];
    
    % Track sequelae
    % 1- schisto light infection
    % 2- schisto heavy infection
    % 3- anemia, low
    % 4- anemia, moderate
    % 5- anemia, heavy
    % 5- schisto mortality
    
    total_DALY_sequelae_a(i,:)= [schisto_prob1_a*schisto_weight(1)* total_pop_a, schisto_prob2_a*schisto_weight(3)* total_pop_a, prob_anemia_low_a* anemia_weight(1)* total_pop_a,  prob_anemia_mod_a* anemia_weight(2)* total_pop_a, prob_anemia_severe_a* anemia_weight(3)* total_pop_a, DALY_matrix_count_a(10)]*(1/12);
    total_DALY_sequelae_c(i,:)= [schisto_prob1_c*schisto_weight(1)* total_pop_c, schisto_prob2_c*schisto_weight(3)* total_pop_c, prob_anemia_low_c* anemia_weight(1)* total_pop_c,  prob_anemia_mod_c* anemia_weight(2)* total_pop_c, prob_anemia_severe_c* anemia_weight(3)* total_pop_c, DALY_matrix_count_c(10)]*(1/12);
    total_DALY_sequelae_pre(i,:)= [schisto_prob1_pre*schisto_weight(1)* total_pop_pre, schisto_prob2_pre*schisto_weight(3)* total_pop_pre, prob_anemia_low_pre* anemia_weight(1)* total_pop_pre,  prob_anemia_mod_pre* anemia_weight(2)* total_pop_pre, prob_anemia_severe_pre* anemia_weight(3)* total_pop_pre, DALY_matrix_count_pre(10)]*(1/12);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Costing analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) MDA_strategy (scalar)
%        (0)None
%        (1)MDA (SAC)
%        (2)MDA (SAC, pre-SAC)
%        (3)MDA (Community)
% 2) Snail_strategy (scalar)
%        (0)None
%        (1)Molluscicides
% 3) MDA_coverage (g_MDA; 1x3 vector; fraction e.g., 75% is "0.75")
%        Note: [pre-SAC, SAC, adults]
% 4) Snail_coverage (g_snail; scalar; fraction e.g., 75% is "0.75") )
% 5) MDA frequency (MDA_freq)
%        (1)Annual
%        (2)Biannual (2/year)
% 6) Molluscicides frequency (snail_freq)
%        (1)Annual
%        (2)Biannual (2/year)
%        (12)Monthly

% MDA (annual cost)
if MDA_strategy==0
    total_cost=0;
elseif MDA_strategy==1
    total_cost=(cost_children*total_pop_c*g_MDA_c) * MDA_freq;
elseif MDA_strategy==3
    total_cost=((cost_pre*total_pop_pre*g_MDA_pre) + (cost_children*total_pop_c*g_MDA_c) + (cost_adults*total_pop_a*g_MDA_a)) * MDA_freq;
end

% Snail control (annual cost)
% Remember to add MDA cost
if snail_strategy==0
    total_cost;
elseif snail_strategy==1
    total_cost= total_cost + cost_snail_control_community*g_snail*snail_freq;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discounting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Apply discounting to costs and disability at 3% annually as per WHO
% guidelines.

% Discount DALYs (monthly time step)
% 3% annual rate is equivalent to 0.24663% monthly rate
for t=1:(time_model-1)
    if discounting_on==1; discount=1/((1+(0.0024663))^(t-1)); end
    if discounting_on==0; discount=1; end
    
    total_DALY_iteration(t+1,1)=total_DALY_iteration(t+1,1)*discount;
    total_DALY_iteration(t+1,2)=total_DALY_iteration(t+1,2)*discount;
    total_DALY_iteration(t+1,3)=total_DALY_iteration(t+1,3)*discount;
end

% Discount costs (yearly time step)
% 3% annual rate
discount_freq=1;

for t=1:(((time_model-1)/12))
    if discounting_on==1; discount=1/((1+disc_rate)^(t-1)); end
    if discounting_on==0; discount=1; end
    
    if t==1 || round((t-1)/discount_freq)==((t-1)/discount_freq)
        
        cost_vector(t)= total_cost*discount;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cost-effectiveness analysis: ICER and final calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Total DALY summary for intervention
total_DALYs=sum(total_DALY_iteration);
total_DALYs=sum(total_DALYs);

% Total DALY summary for intervention by adults and children, respectively
total_DALYs_age=sum(total_DALY_iteration);

total_cost=sum(cost_vector);
total_DALY_sequelae_a=sum(total_DALY_sequelae_a);
total_DALY_sequelae_c=sum(total_DALY_sequelae_c);
total_DALY_sequelae_pre=sum(total_DALY_sequelae_pre);

% Provide total costs and DALYs for ICER calculations
ICER_values=[total_cost, total_DALYs];

total_DALY_iteration=total_DALY_iteration;

end
