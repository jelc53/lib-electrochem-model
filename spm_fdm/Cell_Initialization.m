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

%% SPECIFY: Input current profile
%To interface with the main code, any current profile should be
%defined as a vector titled 'param.I_data'. The corresponding time vector 
%should be defined as 'param.t_data'. 
% Variables 'param.I_data' and 'param.t_data' should be the same lengths. 
% Also, specify 'param.dt' for your current
% profile, which is the time-step size for the 'param.t_data'.

%Current profile data can be loaded from .xlsx or .mat files, as long as
%the current and time points are ultimately defined using the variable
%names outlined above.


%% Specify Cell Cycling

% For a custom input current profile (constant current discharge)
param.t_duration = 3600; % simulation time in seconds
param.dt = 0.5; % sampling time            
param.t_data = [0:param.dt:param.t_duration]'; % time vector
I_multiply = 2; %current magnitude
param.I_data = I_multiply*ones(length(param.t_data),1); % current vector

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

% Average concentration plots
figure()
subplot(2,1,1)
plot(param.t_data/3600, mean(cs_p), 'b','Linewidth',2); grid on
xlabel('Time [h]','Fontsize',16,'interpreter','latex')
ylabel('$c_{avg,p}$ [$mol/m^3$]','Fontsize',16,'interpreter','latex')
title('Average Solid Concentration - Positive Electrode','Fontsize',16,'interpreter','latex')
set(gca,'Fontsize',16)
xlim([0 max(param.t_data/3600)])

subplot(2,1,2)
plot(param.t_data/3600, mean(cs_n), '--r','Linewidth',2); grid on
xlabel('Time [h]','Fontsize',16,'interpreter','latex')
ylabel('$c_{avg,p}$ [$mol/m^3$]','Fontsize',16,'interpreter','latex')
title('Average Solid Concentration - Negative Electrode','Fontsize',16,'interpreter','latex')
set(gca,'Fontsize',16)
xlim([0 max(param.t_data/3600)])
