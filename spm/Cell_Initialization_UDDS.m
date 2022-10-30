%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Simulate a lithium-ion cell cell %%%%%%%%

% % Model description: 
% Electrochemical model: Single Particle Model (ESPM)

% % Specify the following inputs:
% 1) simulation time: t_duration
% 2) current input: I_data
% 3) Number of radial discretization grids in SPPM: Nr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
warning off
clc

%% UDDS Profile

% Note 1: Save the excel file "NMC_Cell_H1_T23_UDDS.xlsx" in your working
% folder

% Note 2: Copy this UDDS profile code section and paste it in your
% SPM Cell_Initialization.m script. Comment out the following variables of
% the code "param.t_duration; param.dt; param.t_data; I_multiply = 2; 
% param.I_data", since we will be replacing them with the UDDS
% counterparts.

% Note 3: The UDDS profile is a charge-depleting profile for an Electric 
% Vehicle. It traverses only around 20% of SOC (or depth of discharge). For
% instane, if we start from 100% SOC, the final SOC at the end of the cycle
% will be around 80%.

data = xlsread('NMC_Cell_H1_T23_UDDS.xlsx');
t_exp = data(:,2); % in seconds
I_exp = -data(:,3); % in Amperes (negative is charging, and positive is discharging)

param.dt = 1; % sampling time   
param.I_data = I_exp;
param.t_data = t_exp;

% % Note that one cycle of UDDS will only discharge the cell by around
% 15-20% of SOC. If you want to traverse through a larger range of SOC
% then concatenate multiple UDDS profiles.

% UDDS_concatenate = 3; % Keep 1 if one UDDS cycle only; max limit = 4;
% param.I_data = [];
% param.t_data = [];
% for i = 1:UDDS_concatenate
% param.I_data = [param.I_data; I_exp];
% param.t_data = [param.t_data; (i-1)*t_exp(end)+t_exp];
% end

% %Interpolate experimental profile to desired timestep size
% t_data = [t_exp(1):param.dt:t_exp(end)]'; % time vector
% % t_data = [t_exp(1):param.dt:180]'; % time vector - shortened
% I_multiply = 2; % Control severity of current profile 
% I_data = I_multiply*interp1(t_exp, I_exp, t_data);
% 
% param.I_data = [];
% param.t_data = [];
% for i = 1:cycles
% param.I_data = [param.I_data; I_data];
% param.t_data = [param.t_data; (i-1)*t_data(end)+t_data];
% end

%% Specify Cell Cycling

% % For a custom input current profile (constant current discharge)
% param.t_duration = 3600; % simulation time in seconds
% param.dt = 1; % sampling time            
% param.t_data = [0:param.dt:param.t_duration]'; % time vector
% I_multiply = 2; %current magnitude
% param.I_data = I_multiply*ones(length(param.t_data),1); % current vector

% Specify number of additional cycles beyond initial charge / discharge
% 'Run_module' script will concatenate alternating charge / discharge
% current profiles as needed to meet the input # of cycles
param.cycles = 0;

% SPECIFY: Finite Difference / Volume Discretization Parameters

param.Nc = 1;      % Number of cells in series (Keep this 1)
param.Nr = 10;     % Number of radial discretization grids in ESPM

% SPECIFY: degree of variation in model parameters and initialize

param.SOC_ref = 1;                            % Initial SOC

% Initialize all model parameters
run ModelParameters.m

% SPECIFY: initial conditions of all cells
% Specify Ambient Temperature
param.T_amb = 298; %[K]  
param.T_cell = param.T_amb;

% Initial lithium concentration in solid phase for all cells corresponding
% to the cell SOC
[cs_initial, csn0, csp0] = conc_initial_sd(param.SOC_ref, param);

% Group all initial state variables in the state vector
x_initial = [cs_initial];

%% RUN SIMULATION
tic

% Call main 'Module_Run' script
[V_cell, soc_bulk_n, soc_bulk_p, cs_n, cs_p, Up, Un, V_oc, param] = Cell_Run(x_initial,param);

run_time = toc

%% The following sections contain code for frequently used plots

% Input current profile plot
figure()
plot(param.t_data/3600, param.I_data, 'b','Linewidth',2); grid on
xlabel('Time [h]','Fontsize',16,'interpreter','latex')
ylabel('Current [A]','Fontsize',16,'interpreter','latex')
title('Current Profile','Fontsize',16,'interpreter','latex')
set(gca,'Fontsize',16)
xlim([0 max(param.t_data/3600)])

% Cell voltage and electrode OCP plots
figure()
subplot(3,1,1)
plot(param.t_data/3600, V_cell, 'b','Linewidth',2); grid on
xlabel('Time [h]','Fontsize',16,'interpreter','latex')
ylabel('Voltage [V]','Fontsize',16,'interpreter','latex')
title('Voltage Response','Fontsize',16,'interpreter','latex')
set(gca,'Fontsize',16)
xlim([0 max(param.t_data/3600)])

subplot(3,1,2)
plot(param.t_data/3600, Up, '--r','Linewidth',2); grid on
xlabel('Time [h]','Fontsize',16,'interpreter','latex')
ylabel('Voltage [V]','Fontsize',16,'interpreter','latex')
title('Open Circuit Potential - Positive Electrode','Fontsize',16,'interpreter','latex')
set(gca,'Fontsize',16)
xlim([0 max(param.t_data/3600)])

subplot(3,1,3)
plot(param.t_data/3600, Un, ':g','Linewidth',2); grid on
xlabel('Time [h]','Fontsize',16,'interpreter','latex')
ylabel('Voltage [V]','Fontsize',16,'interpreter','latex')
title('Open Circuit Potential - Negative Electrode','Fontsize',16,'interpreter','latex')
set(gca,'Fontsize',16)
xlim([0 max(param.t_data/3600)])

% Bulks SOC plots
figure()
subplot(2,1,1)
plot(param.t_data/3600, soc_bulk_p, 'b','Linewidth',2); grid on
xlabel('Time [h]','Fontsize',16,'interpreter','latex')
ylabel('SOC [-]','Fontsize',16,'interpreter','latex')
title('Bulk SOC - Positive Electrode','Fontsize',16,'interpreter','latex')
set(gca,'Fontsize',16)
xlim([0 max(param.t_data/3600)])

subplot(2,1,2)
plot(param.t_data/3600, soc_bulk_n, '--r','Linewidth',2); grid on
xlabel('Time [h]','Fontsize',16,'interpreter','latex')
ylabel('SOC [-]','Fontsize',16,'interpreter','latex')
title('Bulk SOC - Negative Electrode','Fontsize',16,'interpreter','latex')
set(gca,'Fontsize',16)
xlim([0 max(param.t_data/3600)])

% Surface concentration plots
figure()
subplot(2,1,1)
plot(param.t_data/3600, cs_p(end,:), 'b','Linewidth',2); grid on
xlabel('Time [h]','Fontsize',16,'interpreter','latex')
ylabel('$c_{surf,p}$ [$mol/m^3$]','Fontsize',16,'interpreter','latex')
title('Surface Concentration - Positive Electrode','Fontsize',16,'interpreter','latex')
set(gca,'Fontsize',16)
xlim([0 max(param.t_data/3600)])

subplot(2,1,2)
plot(param.t_data/3600, cs_n(end,:), '--r','Linewidth',2); grid on
xlabel('Time [h]','Fontsize',16,'interpreter','latex')
ylabel('$c_{surf,n}$ [$mol/m^3$]','Fontsize',16,'interpreter','latex')
title('Surface Concentration - Negative Electrode','Fontsize',16,'interpreter','latex')
set(gca,'Fontsize',16)
xlim([0 max(param.t_data/3600)])
