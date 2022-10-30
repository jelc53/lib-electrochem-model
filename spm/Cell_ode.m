function [dx_dt] = Cell_ode(t_in, x_in, param)
% Solving solid phase PDE using ODE solver
% Interpolate current
input_crt = interp1(param.t_data, param.I_data, t_in);   

%% Separate electrochemical, thermal & aging state variables from x_in vector
%Define solid concentrations
cs = x_in(1:(param.Nr-1)*param.Nc*2);           %Solid Concentration states
cs_n = cs(1:(param.Nr-1)*param.Nc);             %Anode Solid Conc.
cs_p = cs((param.Nr-1)*param.Nc+1:end);         %Cathode Solid Conc. 

%% Surface concentration -> surface stoichiometry

%Allocate memory
theta_surf_n(param.Nc,1) = 0;
theta_surf_p(param.Nc,1) = 0;

for i = 1:param.Nc
    % Surface concentration -> surface stoichiometry
    theta_surf_n(i,:) = cs_n(i*(param.Nr-1))/param.c_n_max;
    theta_surf_p(i,:) = cs_p(i*(param.Nr-1))/param.c_p_max;  
end


%% Open circuit potential and overpotential
ocp_p = U_p(theta_surf_p);
ocp_n = U_n(theta_surf_n);
eta_p = eta_cathode(theta_surf_p, param.ce0, param.T_cell, input_crt, param, param.kp_ref);
eta_n = eta_anode(theta_surf_n, param.ce0, param.T_cell, input_crt, param, param.kn_ref);

% Cell voltages
V_cell = ocp_p - ocp_n + eta_p - eta_n - param.R_l*input_crt;

% Heat generation in each cell
V_loss = (ocp_p - ocp_n - V_cell); 

%% Solve solid phase ODEs

% Coefficients of discretized ODEs
alpha_n = param.Dsn_ref/(param.delta_xn^2);
alpha_p = param.Dsp_ref/(param.delta_xp^2);

A_csn = param.A_sd.*alpha_n; %Full Anode Solid Conc. State A-matrix
A_csp = param.A_sd.*alpha_p; %%Full Cathode Solid Conc. State A-matrix

dcsn_dt = A_csn*cs_n + param.B_sd_n.*input_crt;
dcsp_dt = A_csp*cs_p + param.B_sd_p*input_crt;

dcs_dt = [dcsn_dt; dcsp_dt];


%% Final ODE formulation
dx_dt = [dcs_dt];
end

