function J = ESPM_fit(x)
    warning off
    
    global pso_iter err_theta_surf_p err_Up t_data dt V_data I_data SOC_cc flag_up...
           SOC_IC Q_IC Lsei_IC T_amb n_iter x_ rms_values
     
    n_iter = n_iter + 1;
    x_(n_iter,:) = x;

try
    %% Specify Cell Cycling
    % For a custom input current profile
    param.t_duration = t_data(end);                % simulation time in seconds
    param.dt = dt;                                 % sampling time            
    param.t_data = [0:param.dt:param.t_duration]'; % time vector
    param.I_data = I_data;                         % current vector

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
    param.SOC_ref = SOC_IC;                        % Average SOC
    param.SOC_deviation = 0;                       % Average SOC deviation
    param.R_l_deviation = 0;                       % Lumped resistance Stochastic Deviation
    param.Lsei_ref = Lsei_IC;                      % SEI Thickness Reference Value
    param.SEI_deviation = 0;                       % SEI Thickness Stochastic Deviation

    % Initialize all model parameters
    run ModelParameters.m

    % SPECIFY: initial conditions of all cells
    % Specify Ambient Temperature
    param.T_amb = 273 + T_amb; %[K]

    % Induce random (Gaussian) variation around the average SOC
    SOC_cells = soc_variation(param);  

    % Initial lithium concentration in solid phase for all cells corresponding
    % to the cell SOC
    [cs_initial,csn0,csp0] = conc_initial_sd(SOC_cells, param);

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

    %% RUN SIMULATION
    [V_cell, T_core, T_surf,soc_bulk_n, soc_bulk_p, cs_n, cs_p,...
     ce_all,L_sei,i_s,Q,V_oc,R_el,R_sei,Csolv,aina_n,aina_p,...
     c_sei,c_li,i_lpl,L_film,af_n,af_p,param] = Module_Run(x_initial,param);
 
    % Fit
    V_cell_fit = interp1(param.t_data,V_cell,t_data,'linear','extrap');
    SOCp_fit = interp1(param.t_data,soc_bulk_p,t_data,'linear','extrap');
    SOCn_fit = interp1(param.t_data,soc_bulk_n,t_data,'linear','extrap');

    %% Compute cost function 
    if sum(isreal(V_cell)) == 0 
       J = 10^2;
    else
       % Cost function
       J_v = rms(V_cell_fit-V_data);
       J_SOCp = rms(SOC_cc-SOCp_fit');
       J_SOCn = rms(SOC_cc-SOCn_fit');
       J = J_v + J_SOCp + J_SOCn;
    end
catch
    J = 10^2;
end

% Save
rms_values(n_iter) = J;
save(['opt_res/hist/RMS_opt_' num2str(pso_iter) '.mat'],'rms_values','x_')
end