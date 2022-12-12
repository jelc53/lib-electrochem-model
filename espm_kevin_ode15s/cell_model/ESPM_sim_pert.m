function [V_cell, T_core, T_surf,soc_bulk_n, soc_bulk_p, cs_n, cs_p,...
          ce_all, L_sei,i_s,Q,V_oc,R_el,R_sei,Csolv,aina_n,aina_p,c_sei,c_li,i_lpl,L_film,af_n,af_p,param] = ESPM_sim_pert(x,dt,t_data,I_data,SOC_IC,Q_IC,Lsei_IC,T_amb,pct,theta_names,theta_mask) 
    
    % Function to produce sensitivity analysis on outputs by
    % varying parameters according to input mask and percentage
    % perturbation around the previously-identified parameters.
    
    % Function inputs:
    % x == vector of previously-identified parameters
    % dt == sampling time of experimental data
    % t_data == time vector from experimental data
    % I_data == current vector from experimental data (nominally Amperes)
    % SOC_IC == State of Charge Initial Condition (range [0, 1])
    % Q_IC == Cell Capacity Initial Condition (nominally Ampere-hours)
    % Lsei_IC == SEI Layer Thickness Initial Condition (nominally meters?)
    % T_amb == Ambient temperature of experimental data
    % pct == Percentage variation to induce in all/each parameters (range [0.1])
    % theta_names == Cell vector of parameter names to perturb
    % theta_mask == Mask of parameters to vary with this function call (length == length(theta_names))
    
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
    param.SOC_ref = SOC_IC;                         % Average SOC
    param.SOC_deviation = 0;                        % Average SOC deviation
    param.R_l_deviation = 0;                        % Lumped resistance Stochastic Deviation
    param.Lsei_ref = Lsei_IC;                       % SEI Thickness Reference Value
    param.SEI_deviation = 0;                        % SEI Thickness Stochastic Deviation

    % Initialize all model parameters
    run ModelParameters_geom.m
    
    % ------ NEW CODE ------

    % Apply percentage pct on mask theta_mask to subset of nominal 
    % parameters named in theta_names.

    % Mask indicates whether or not we are applying the perturbation:
    % 1 indicates we add pct to the value of the parameter 
    % ---> perturbed_param_value = (nominal_parameter_value * (1+pct);
    % 0 indicates that we keep the parameter value as nominal

    for i = 1:length(theta_names)
        
%         % Debugging: Display original value of theta:
%         disp([theta_names{i}, ' Original value: ', num2str(param.(theta_names{i}))])

        % Compute pct deviation for theta;
        delta_theta = theta_mask(i)*pct;

        if strcmp(theta_names{i}, 'theta100_p') % Avoid undervoltage issue with theta100_p
            delta_theta = -delta_theta; 
        end

        % Set new perturbed theta value:
%         param = setfield(param, theta_names{i}, getfield(param, theta_names{i}) * (1 + delta_theta));
        param.(theta_names{i}) = param.(theta_names{i}) * (1 + delta_theta);

%         % Debugging: Display new perturbed value of theta:
%         disp([theta_names{i}, ' Perturbed value: ', num2str(param.(theta_names{i}))])
    end

    % Recompute dependent parameters (e.g. lumped resistance, electrolyte porosity, surface area, grid interval sizes)

    param.R_l = param.R_l_ref*ones(param.Nc,1) + (param.R_l_deviation)*randn(param.Nc,1);

    % Electrolyte porosity
    param.eps_filler_n = 0.75-param.epsilon_n;
    param.eps_filler_p = 0.665-param.epsilon_p; 
    param.eps_el_n_ref = 1-param.epsilon_n-param.eps_filler_n;
    param.eps_el_p = 1-param.epsilon_p-param.eps_filler_p;
    param.eps_el_s = 0.47; 
    
    % Specific interfacial (electroactive) surface area for anode and cathode
    param.a_sn = 3*param.epsilon_n/param.Rs_n;   
    param.a_sp = 3*param.epsilon_p/param.Rs_p; 
    
    % For finite volume method
    param.delta_n = param.Ln/(param.Nx_n);
    param.delta_s = param.Ls/(param.Nx_s);
    param.delta_p = param.Lp/(param.Nx_p);

    % Derived parameters for solid phase discretization
    % Grid interval size for anode and cathode particle
    param.delta_xn = param.Rs_n/(param.Nr-1);    
    param.delta_xp = param.Rs_p/(param.Nr-1);  

    % ------ END OF NEW CODE ------



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
           
    %% RUN SIMULATION
    [V_cell, T_core, T_surf,soc_bulk_n, soc_bulk_p, cs_n, cs_p,...
        ce_all, L_sei,i_s,Q,V_oc,R_el,R_sei,Csolv,aina_n,aina_p,c_sei,c_li,i_lpl,L_film,af_n,af_p,param] = Module_Run(x_initial,param);

end