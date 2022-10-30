%% UDDS Profile

% Note 1: Save the excel file "NMC_Cell_H1_T23_UDDS.xlsx" in your working
% folder.

% Note 2: Copy this UDDS profile code section and paste it in your
% SPM Cell_Initialization.m script. Comment out the following variables of
% the old code "param.t_duration; param.dt; param.t_data; I_multiply; 
% param.I_data" since we will be replacing them with their UDDS
% counterparts below.

% Note 3: The UDDS profile is a charge-depleting profile for an Electric 
% Vehicle. It traverses only around 20% of SOC (or depth of discharge). For
% instance, if we start from 100% SOC, the final SOC at the end of the cycle
% will be around 80%.

data = xlsread('NMC_Cell_H1_T23_UDDS.xlsx');
t_exp = data(:,2); % in seconds
I_exp = -data(:,3); % in Amperes (negative is charging, and positive is discharging)

param.dt = 0.1; % sampling time   
param.I_data = I_exp;
param.t_data = t_exp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%