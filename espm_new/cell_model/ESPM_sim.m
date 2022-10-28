function [V_cell, R_l, T_core, T_surf,soc_bulk_n, soc_bulk_p, cs_n, cs_p,...
          ce_all, L_sei,i_s,Q,V_oc,R_el,R_sei,Csolv,aina_n,aina_p,...
          c_sei,c_li,i_lpl,L_film,af_n,af_p,param] = ESPM_sim(x, dt, t_data, I_data, ...
          SOC_cc, SOC_IC,Q_IC,Lsei_IC,T_amb) 

    % For a custom input current profile
%     param.t_duration = t_data(end);                % simulation time in seconds
%     param.dt = dt;                                 % sampling time            
%     param.t_data = [0:param.dt:param.t_duration]'; % time vector
    param.t_data = t_data;
    param.I_data = I_data;                         % current vector
    param.SOC_cc = SOC_cc;                         % Coulomb-counted SOC (needed for lumped resistance)

    % Specify number of additional cycles beyond initial charge / discharge
    % 'Run_module' script will concatenate alternating charge / discharge
    % current profiles as needed to meet the input # of cycles
    param.cycles = 0;

    % SPECIFY: Finite Difference / Volume Discretization Parameters
    param.Nc = 1;       % Number of cells in series
    param.Nr = 10;      % Number of radial discretization grids in ESPM
    param.Nx_n = 10;    % Number of cartesian discretization grids in ESPM
    param.Nx_s = 10;    % Number of cartesian discretization grids in ESPM
    param.Nx_p = 10;    % Number of cartesian discretization grids in ESPM
    param.Nsei = 10;    % SEI Layer Discretization

    % SPECIFY: degree of variation in model parameters and initialize
    param.SOC_ref = SOC_IC;                         % Average SOC
    
    param.SOC_deviation = 0;                        % Average SOC deviation
    param.R_l_deviation = 0;                        % Lumped resistance Stochastic Deviation
    param.Lsei_ref = Lsei_IC;                       % SEI Thickness Reference Value
    param.SEI_deviation = 0;                        % SEI Thickness Stochastic Deviation

    % Initialize all model parameters
    run ModelParameters.m

    % SPECIFY: initial conditions of all cells
    % Specify Ambient Temperature
    param.T_amb = 273 + T_amb; %[K]

    % Induce random (Gaussian) variation around the average SOC
    SOC_cells = soc_variation(param);  

    % Initial lithium concentration in solid phase for all cells corresponding
    % to the cell SOC
    [cs_initial,csn0,csp0] = conc_initial_sd(SOC_cells,param);

    % Initial lithium concentration in electrolyte phase for all cells
    param.ce_states = param.Nx_n + param.Nx_s + param.Nx_p;
    ce_initial = param.ce0*ones(param.Nc*param.ce_states,1);

    % Initial temperature for all cells
    param.T_states = 2*param.Nc; % Factor of 2 because of two-state model
    T_initial = param.T_amb*ones(param.T_states,1); % Same temperature as ambient

    % Initial SEI layer thickness for all cells (L_sei)
    sei_initial = sei_variation(param);

    % Initial Capacity for all cells
    Q_initial = Q_IC*ones(param.Nc,1);
    
    % Initial inactive surface area (NB: no variation w.r.t. to cartesian/radial coordinates)
    % Two states: one for anode, one for cathode
    aina_initial_n = param.aina_ic_n*ones(param.Nc,1);
    aina_initial_p = param.aina_ic_p*ones(param.Nc,1);
    
    % Initial inactive surface area (NB: no variation w.r.t. to cartesian/radial coordinates)
    % Two states: one for anode, one for cathode
    c_sei_initial = zeros(param.Nc,1);
    c_li_initial = zeros(param.Nc,1);

    % Initial Solvent Concentration in SEI Layer
    Csolv_initial = param.eps_sei*param.c_solv*ones(param.Nc*param.Nsei,1);
    
    % Initial film Layer
    film_initial = sei_initial;
 
    % Group all initial state variables
    x_initial = [cs_initial; ce_initial; T_initial; sei_initial; Q_initial; aina_initial_n; aina_initial_p; c_sei_initial; c_li_initial; film_initial; Csolv_initial];
    x_initial = [cs_initial;ce_initial];
                 
    %% RUN SIMULATION
 
    [V_cell, R_l, T_core, T_surf,soc_bulk_n, soc_bulk_p, cs_n, cs_p,...
        ce_all, L_sei,i_s,Q,V_oc,R_el,R_sei,Csolv,aina_n,aina_p,c_sei,c_li,i_lpl,L_film,af_n,af_p,param] = Module_Run(x_initial,param);

end