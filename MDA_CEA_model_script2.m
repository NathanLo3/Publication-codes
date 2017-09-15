%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to run "MDA_CEA_model_CaseWestern.m"
% Written by Nathan Lo
%
% Last updated: 8/16/17
%
% Collaboration with Case Western
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
% 7) excel_decision
%        (1) Which sheet (1, 2, 3...)
%        (2) Base case (0), lower limit (-1), upper limit (1)

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

% Case Western data inputs (based on tested strategies)
MDA_strategy_input=[0, 3, 3, 1, 1, 0, 0, 3, 3, 3, 3, 1, 1, 1, 1];
snail_strategy_input=[0, 0, 0, 0, 0, 1, 1, 1, 1 ,1 ,1 ,1 ,1 ,1 ,1];
g_MDA_input=[0.75 0.75 0.75];
g_snail_input=1;
MDA_freq_input=[0, 1, 2, 1, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 2];
snail_freq_input=[0, 0, 0, 0, 0, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2];

for j=1:3
    which_analysis=[0 -1 1];
    % 0- base, -1- lower bounds, 1- upper bounds
    excel_decision(2)=which_analysis(j);
    
    for i=1:length(MDA_strategy_input)
        
        excel_decision(1)=i;
        MDA_strategy= MDA_strategy_input(i);
        snail_strategy= snail_strategy_input(i);
        g_MDA= g_MDA_input;
        g_snail= g_snail_input;
        MDA_freq= MDA_freq_input(i);
        snail_freq= snail_freq_input(i);
        
        [total_DALY_sequelae_pre, total_DALY_sequelae_c, total_DALY_sequelae_a, ICER_values, total_DALYs_age, total_cost, total_DALY_iteration] = MDA_CEA_model_CaseWestern(MDA_strategy, snail_strategy, g_MDA, g_snail, MDA_freq, snail_freq, excel_decision);
        
        if excel_decision(2)==0
            sim_data_base(i,:)=[i-1, ICER_values];
        elseif excel_decision(2)==-1
            sim_data_lower(i,:)=[i-1, ICER_values];
        elseif excel_decision(2)==1
            sim_data_upper(i,:)=[i-1, ICER_values];
        end
        
    end
end

sim_data_combined=[sim_data_base sim_data_lower sim_data_upper];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute cost-effectiveness (ICER)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Columns
% 1- Strategy number
% 2- Cost (US$)
% 3- Disbility (DALYs)

CEA_data_sorted_base=sortrows(sim_data_base,-3);
CEA_data_sorted_base_temp=sortrows(sim_data_base,-3);
CEA_data_sorted_lower=sortrows(sim_data_lower,-3);
CEA_data_sorted_lower_temp=sortrows(sim_data_lower,-3);
CEA_data_sorted_upper=sortrows(sim_data_upper,-3);
CEA_data_sorted_upper_temp=sortrows(sim_data_upper,-3);

% Create Dominated (big D) vector
for i=1:length(CEA_data_sorted_base_temp)
    dominant_vector_base(i)= sum(CEA_data_sorted_base_temp(i,2)>CEA_data_sorted_base_temp(i+1:length(CEA_data_sorted_base_temp),2));
    dominant_vector_lower(i)= sum(CEA_data_sorted_lower_temp(i,2)>CEA_data_sorted_lower_temp(i+1:length(CEA_data_sorted_lower_temp),2));
    dominant_vector_upper(i)= sum(CEA_data_sorted_upper_temp(i,2)>CEA_data_sorted_upper_temp(i+1:length(CEA_data_sorted_upper_temp),2));
end
% element 0 is non-dominated
% element>0 is Dominated

% ICER_dominated = zeros(1,length(CEA_data_sorted_base_temp(:,1)))+1;
% % Compute ICER, removing dominated comparisons
% for i=2:length(CEA_data_sorted_base_temp)
%
%     vec = CEA_data_sorted_base_temp(i,2) < CEA_data_sorted_base_temp(1:(i-1),2);
%
%     if sum(vec)>0
%         ICER_dominated(find(vec==1))=0;
%     end
%
% end
track_base=1;
track_lower=1;
track_upper=1;
for i=2:length(CEA_data_sorted_base_temp)
    
    % Base
    if dominant_vector_base(i)==0
        
        while(track_base<i)
            rel_iter=max(find(dominant_vector_base(1:i-1)==0));
            rel_iter;
            ICER_calc_base(i)= [CEA_data_sorted_base_temp(i,2)-CEA_data_sorted_base_temp(rel_iter,2)] / [CEA_data_sorted_base_temp(rel_iter,3) - CEA_data_sorted_base_temp(i,3)];
            
            if sum(ICER_calc_base(i)<ICER_calc_base(1:i-1))>0
                iter_change = max(find(ICER_calc_base(i)<ICER_calc_base(1:i-1)));
                dominant_vector_base(iter_change)=1;
                ICER_calc_base(iter_change)=0;
            else
                track_base= track_base+1;
            end
        end
    else
        track_base= track_base+1;
        ICER_calc_base(i)=0;
    end
    
    
    
    % lower
    if dominant_vector_lower(i)==0
        
        while(track_lower<i)
            rel_iter=max(find(dominant_vector_lower(1:i-1)==0));
            rel_iter;
            ICER_calc_lower(i)= [CEA_data_sorted_lower_temp(i,2)-CEA_data_sorted_lower_temp(rel_iter,2)] / [CEA_data_sorted_lower_temp(rel_iter,3) - CEA_data_sorted_lower_temp(i,3)];
            
            if sum(ICER_calc_lower(i)<ICER_calc_lower(1:i-1))>0
                iter_change = max(find(ICER_calc_lower(i)<ICER_calc_lower(1:i-1)));
                dominant_vector_lower(iter_change)=1;
                ICER_calc_lower(iter_change)=0;
            else
                track_lower= track_lower+1;
            end
        end
    else
        track_lower= track_lower+1;
        ICER_calc_lower(i)=0;
    end
    
    
    
    % upper
    if dominant_vector_upper(i)==0
        
        while(track_upper<i)
            rel_iter=max(find(dominant_vector_upper(1:i-1)==0));
            rel_iter;
            ICER_calc_upper(i)= [CEA_data_sorted_upper_temp(i,2)-CEA_data_sorted_upper_temp(rel_iter,2)] / [CEA_data_sorted_upper_temp(rel_iter,3) - CEA_data_sorted_upper_temp(i,3)];
            
            if sum(ICER_calc_upper(i)<ICER_calc_upper(1:i-1))>0
                iter_change = max(find(ICER_calc_upper(i)<ICER_calc_upper(1:i-1)));
                dominant_vector_upper(iter_change)=1;
                ICER_calc_upper(iter_change)=0;
            else
                track_upper= track_upper+1;
            end
        end
    else
        track_upper= track_upper+1;
        ICER_calc_upper(i)=0;
    end
end






CEA_data_final= [CEA_data_sorted_base ICER_calc_base' CEA_data_sorted_lower ICER_calc_lower' CEA_data_sorted_upper ICER_calc_upper']

%header_1= {'base case', 'base case', 'base case', 'base case', 'lower bounds', 'lower bounds', 'lower bounds', 'lower bounds', 'upper bounds', 'upper bounds', 'upper bounds', 'upper bounds'};
%header_2= {'Strategy', 'Total cost (US$)', 'Total disability (DALYs)', 'ICER (US$/DALY)', 'Strategy', 'Total cost (US$)', 'Total disability (DALYs)', 'ICER (US$/DALY)', 'Strategy', 'Total cost (US$)', 'Total disability (DALYs)', 'ICER (US$/DALY)'};
%A= [header_1; header_2];
B= CEA_data_final;
filename = 'Base_case_090217_lowprev.xls';
%filename = 'Scenario_090917_highprev_8.xls';
xlswrite(filename,B)

