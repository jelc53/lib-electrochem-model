%% Initialize Model Parameters (SI Units)
% Framework for geometric parameter identification
% FOR ELECTROCHEMICAL & THERMAL MODELS: 
% Cell Name: INR21700
% Manufacturer: LG Chem
% Rated Capacity = 4.85 A-h
% Nominal Voltage = 3.63V 
% Cathode Chemistry = NMC / Anode Chemistry = Graphite

param.Rg = 8.314;                   % Gas constant
param.F = 96487;                    % Faraday's constant
param.alpha_cell = 0.5;             % Anode/Cathode transfer coefficient
param.ce0 = 1000;                   % Average electrolyte concentration [mol/m^3]
param.capacity = 4.85;              % Nominal cell capacity (Ah)   

param.c_n_max = 29583;              % Maximum anode concentration (from Chen et al. 2020)
param.c_p_max = 51765;              % Maximum cathode concentration (from Chen et al. 2020)

%% GEOMETRIC PARAMETERS
param.Rs_n = 10^(-x(1));            % Anode particle radius
param.Rs_p = 10^(-x(2));            % Cathode particle radius         
param.A = (0.1027 + 0.10465)/2;     % Area (averaged from anode/cathode area from Chen et al. 2020)
param.Ln = 85.2e-06;                % Anode thickness (from Chen et al. 2020)
param.Lp = 75.6e-06;                % Cathode thickness (from Chen et al. 2020)
param.epsilon_n = x(3);             % Anode solid phase volume fraction (1-poro)
param.epsilon_p = x(4);             % Cathode solid phase volume fraction (1-poro)
param.theta100_n = x(5);            % Anode stoichiometry at 100% SOC
param.theta100_p = x(6);            % Cathode stoichiometry at 100% SOC
param.theta0_n = x(7);              % Anode stoichiometry at 0% SOC
param.theta0_p = x(8);              % Cathode stoichiometry at 0% SOC

% Electrolyte porosity
param.eps_filler_n = 0.75-param.epsilon_n;          % (from Chen et al. 2020)
param.eps_filler_p = 0.665-param.epsilon_p;         % (from Chen et al. 2020)
param.eps_el_n_ref = 1-param.epsilon_n-param.eps_filler_n;
param.eps_el_p = 1-param.epsilon_p-param.eps_filler_p;
param.eps_el_s = 0.47;                              % (from Chen et al. 2020)

% Specific interfacial (electroactive) surface area for anode and cathode
param.a_sn = 3*param.epsilon_n/param.Rs_n;   
param.a_sp = 3*param.epsilon_p/param.Rs_p; 

% Derived parameters for solid phase discretization
% Grid interval size for anode and cathode particle
param.delta_xn = param.Rs_n/(param.Nr-1);    
param.delta_xp = param.Rs_p/(param.Nr-1);  

%% TRANSPORT/KINETIC PARAMETERS
% Diffusion coefficients
param.Dsn_ref = 10^(-x(9));   % Anode diffusion coefficient
param.Dsp_ref = 10^(-x(10));   % Cathode diffusion coefficient

% Reaction rates
param.kn_ref = 10^(-x(11))/param.F;     % Anode reaction rate constant 
param.kp_ref = 10^(-x(12))/param.F;     % Cathode reaction rate constant

% Lumped Contact Resistance
param.R_l = 0.03;                       % Dummy value; not used in this model at this time

% Additional Electrolyte Parameters
param.t0 = 0.2594;                   % Transference Number (from Chen et al. 2020)
param.Ls = 12e-6;                    % Separator Thickness [m] (from Chen et al. 2020)

% Bruggeman coefficients (from Chen et al. 2020) and supported by others (e.g. Prada)
param.brugg_n = 1.5; 
param.brugg_s = 1.5;
param.brugg_p = 1.5;

% For finite volume method
param.delta_n = param.Ln/(param.Nx_n);
param.delta_s = param.Ls/(param.Nx_s);
param.delta_p = param.Lp/(param.Nx_p);


%% THERMAL MODEL PARAMETERS
param.Tref = 298;                   % Reference Temperature [K]
param.Rc = 1.43;                    % Conduction resistance between core and surface [K/W]
param.Cs = 17.15;                   % Heat Capacity of cell core [J/K]
param.Cc = 59.39;                   % Heat Capacity of cell surface [J/K]
param.Ru = 2.6;                     % Convective resistance between surface and surroundings [K/W]
param.Rm = 19.5;                    % Cell to cell heat transfer resistance in module [K/W]

%% AGING PARAMETERS

% SEI
param.k_f = 0*1.18e-22;             % Solvent Reduction Kinetic Constant
param.M_sei = 0.162;                % SEI Layer Molar Mass
param.rho_sei = 1690;               % SEI Layer Density
param.kappa_sei = 1/1;              % SEI Ionic Conductivity
param.beta_val = 0.5;               % Charge Transfer Coefficient for Side Reaction
param.c_solv = 4541;                % Solvent Bulk Concentration
param.Ea_kf = 6e4;                  % Activation Energy of Side Reaction
param.U_s = 0.4;                    % Solvent Reduction Equilibrium Potential
param.eps_sei = 0.05;               % SEI Layer Porosity
param.Dsolv_ref = 0.25*8.84*10^-20; % Solvent Diffusivity in SEI Layer
param.Ea_Dsolv = 5.55*10^-4;        % Activation Energy of Solvent Diffusion
param.delta_xi = 1/(param.Nsei-1);

% LAM
param.flag_aina_p = 0;                   % Flag LAM cathode
param.flag_aina_n = 0;                   % Flag LAM anode
param.aina_ic_n = 0;                     % Inactive surface area initial condition [m^2/m^3]
param.aina_ic_p = 0;                     % Inactive surface area initial condition [m^2/m^3]
param.af_tinit = 0;                      % Fracture time initial condition [s]
param.beta_lam_n = param.flag_aina_n*0;  % Inactive particle evolution coefficient - anode [-]
param.k_lam_n = param.flag_aina_n*0;     % Fracture evolution coefficient - anode [-]
param.beta_lam_p = param.flag_aina_p*0;  % Inactive particle evolution coefficient -cathode [-]
param.k_lam_p = param.flag_aina_p*0;     % Fracture evolution coefficient - cathode [-]

% Li-plating
param.beta_li2sei = 0;              % Li transformed into SEI [-]
param.M_li = 6.94e-3;               % Li molar mass [kg/mol]
param.rho_li = 534;                 % Li density [kg/m^3]
param.alpha_li = 0.5;               % Li charge transfer coefficient [-]
param.i0_li = 0;                    % Li deposition exchange current density [A/m^2]