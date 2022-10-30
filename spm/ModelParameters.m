%% Initialize Model Parameters (SI Units)
% Defining a structure named "params" that contains all parameter values
% BATTERY DETAILS (i.e. geometry, chemistry)
%
% FOR ELECTROCHEMICAL & THERMAL MODELS: 
% Cell Name: US18650VTC4
% Manufacturer: Sony
% Rated Capacity = 2 A-h
% Nominal Voltage = 3.7V / Maximum Voltage = 4.2V / Minimum Voltage = 2.5V
% Cathode Chemistry = NMC / Anode Chemistry = Graphite

param.Rg = 8.314;                 % Gas constant
param.F = 96487;                  % Faraday's constant
param.alpha_cell = 0.5;           % Anode/Cathode transfer coefficient
param.ce0 = 1200;                 % Average electrolyte concentration [mol/m^3]
param.capacity = 2;               % Nominal cell capacity (Ah)   

param.c_n_max = 31080;            % Maximum anode concentration
param.c_p_max = 51830;            % Maximum cathode concentration

% param.Dsn_ref = 2.87e-14;         % Anode diffusion coefficient
% param.Dsp_ref = 4.85e-14;         % Cathode diffusion coefficient
param.Dsn_ref = 1.4493e-14;         % Anode diffusion coefficient
param.Dsp_ref = 5.8811e-14;         % Cathode diffusion coefficient

param.Rs_n = 4.2642e-06;                % Anode particle radius
param.Rs_p = 5.1592e-07;                % Cathode particle radius
param.A = 0.1048;                  % Area
param.Ln = 4.0000e-05;                 % Anode thickness
param.Lp = 3.6550e-05;                 % Cathode thickness
param.epsilon_n = 0.6187;          % Anode solid phase volume fraction
param.epsilon_p = 0.5800;          % Cathode solid phase volume fraction
param.theta100_n = 0.928;         % Anode stoichiometry at 100% SOC
param.theta100_p = 0.3504;        % Cathode stoichiometry at 100% SOC

% param.kn_ref = 3.48e-10;          % Anode reaction rate constant
% param.kp_ref = 4.164e-10;         % Cathode reaction rate constant

%Test - Anirudh's values
param.kn_ref = 1.0e-10;          % Anode reaction rate constant
param.kp_ref = 1.0e-10;         % Cathode reaction rate constant

param.theta0_n = 0.002;           % Anode stoichiometry at 0% SOC
param.theta0_p = 0.9986;          % Cathode stoichiometry at 0% SOC

param.R_l = 0.032;            % Lumped resistance

% ESPM paramters - TRY USING THIS - TEST
% param.Ln = 40e-6;                   % Anode Thickness [m]
param.Ls = 2.5000e-05;                   % Separator Thickness [m]
% param.Lp = 36.55e-6;                % Cathode Thickness [m]
% param.Rs_n = 4.2642e-6;
% param.Rs_p = 5.1592e-7;
% param.A = 0.1048;
% param.epsilon_n = 0.6187;
% param.epsilon_p = 0.5800;
% param.c_n_max = 31080;            % Maximum anode concentration
% param.c_p_max = 51830;            % Maximum cathode concentration
% ESPM parameters - TRY USING THIS - TEST


% Specific interfacial (electroactive) surface area for anode and cathode
param.a_sn = 3*param.epsilon_n/param.Rs_n;   
param.a_sp = 3*param.epsilon_p/param.Rs_p; 

% % Derived parameters for solid phase discretization
% % Grid interval size for anode and cathode particle
param.delta_xn = param.Rs_n/(param.Nr-1);    
param.delta_xp = param.Rs_p/(param.Nr-1);  


%% Following parameters Not Required : IGNORE

% % Electrolyte porosity
% param.eps_filler_n = 0.038;
% param.eps_filler_p = 0.12; %Value from ESPM code
% % param.eps_filler_p = 0.012; %Original value in module code
% % param.eps_filler_s = 0.4;
% param.eps_el_n_ref = 1-param.epsilon_n-param.eps_filler_n;
% param.eps_el_p = 1-param.epsilon_p-param.eps_filler_p;
% 
% % param.eps_filler_s = 0.4;
% % param.eps_el_s = 1-param.eps_filler_s; %NOTE: VERIFY THIS
% 
% param.eps_el_s = 0.4; %NOTE: Switched this with value above
% 
% %Additional Electrolyte Parameters
% %NOTE: NEED TO VERIFY IF THIS SEPARATOR THICKNESS WAS ALSO IDENTIFIED W/
% %ELECTRODE THICKNESSES SPECIFIED ABOVE
% param.t0 = 0.369;                   % Transference Number
% param.brugg = 1.5;                  % Bruggeman coefficient
% param.Ln_el = 40e-6;                   % Anode Thickness [m]
% param.Ls = 25e-6;                   % Separator Thickness [m]
% param.Lp_el = 36.55e-6;                % Cathode Thickness [m]
% 
% %TEST
% % param.Ln = 40e-6;                   % Anode Thickness [m]
% % param.Ls = 25e-6;                   % Separator Thickness [m]
% % param.Lp = 36.55e-6;                % Cathode Thickness [m]
% 
% param.brugg_n = param.brugg;
% param.brugg_s = param.brugg;
% param.brugg_p = param.brugg;
% 
% % Derived parameters for electrolyte phase discretization
% % Grid interval size for anode, cathode and separator regions
% 
% %For finite difference method:
% % param.delta_n = param.Ln/(param.Nx_n -1);
% % param.delta_s = param.Ls/(param.Nx_s -1);
% % param.delta_p = param.Lp/(param.Nx_p -1);
% 
% %For finite volume method:
% param.delta_n = param.Ln/(param.Nx_n);
% param.delta_s = param.Ls/(param.Nx_s);
% param.delta_p = param.Lp/(param.Nx_p);
% 
% % Activation energy for temperature Arrhenius dependency [J/mol]
% param.Ea_Dsp = 25000;
% param.Ea_Dsn = 50000;
% param.Ea_kp = 30000;
% param.Ea_kn = 30000;
% 
% % Thermal Model parameters
% param.Tref = 298;                   % Reference Temperature [K]
% param.Rc = 1.43;                    % Conduction resistance between core and surface [K/W]
% param.Cs = 17.15;                     % Heat Capacity of cell core [J/K]
% param.Cc = 59.39;                    % Heat Capacity of cell surface [J/K]
% param.Ru = 2.6;                    % Convective resistance between surface and surroundings [K/W]
% param.Rm = 19.5;                    % Cell to cell heat transfer resistance in module [K/W]
% % param.Rm = 1.5; %[K/W] - is this a more accurate value to use?
% 
% 
% % Derived parameters for solid phase discretization
% % Grid interval size for anode and cathode particle
% param.delta_xn = param.Rs_n/(param.Nr-1);    
% param.delta_xp = param.Rs_p/(param.Nr-1);  
% 
% 
% %FOR SENSITIVITY ANALYSIS:
% param.De_scaling = 1; %Electrolyte Diffusivity Scaling Factor
% param.K_el_scaling = 1; %Electrolyte Conductivity Scaling Factor
% 
% 
% % Aging parameters
% param.k_f = 1.18e-22;               % Solvent Reduction Kinetic Constant
% param.M_sei = 0.162;                % SEI Layer Molar Mass
% param.rho_sei = 1690;               % SEI Layer Density
% param.kappa_sei = 17.5e-5;          % SEI Ionic Conductivity
% param.beta_val = 0.5;               % Charge Transfer Coefficient for Side Reaction
% param.c_solv = 4541;                % Solvent Bulk Concentration
% param.Ea_kf = 6e4;                  % Activation Energy of Side Reaction
% param.U_s  = 0.4;                   % Solvent Reduction Equilibrium Potential
% % NEW parameters for Solvent Diffusion Model
% param.eps_sei = 0.05;                % SEI Layer Porosity
% % param.Dsolv_ref = 8.84*10^-20;      % Solvent Diffusivity in SEI Layer
% param.Dsolv_ref = 0.25*8.84*10^-20;      % Solvent Diffusivity in SEI Layer
% param.Ea_Dsolv = 5.55*10^-4;        % Activation Energy of Solvent Diffusion
% param.delta_xi = 1/(param.Nsei-1);