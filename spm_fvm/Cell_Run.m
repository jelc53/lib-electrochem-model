function[V_cell, soc_bulk_n, soc_bulk_p, ...
    cs_n, cs_p, ocp_p, ocp_n, V_oc, param] = Cell_Run(x_initial,param)

%% Generate matrices for solid phase, electrolyte phase, thermal, and aging
% Solid phase discretization matrices
[param.A_sd, param.B_sd_n, param.B_sd_p] = matrices_solidphase(param);

%% Solve ODEs 
tspan = param.t_data;
reltol=5.0e-05; abstol=5.0e-05;
event_formatted = @(t,x) physical_event_function(t,x,param); 
options=odeset('RelTol',reltol,'AbsTol',abstol,'Events',event_formatted);
% [t_out, x_out, te, xe, ie] = ode45(@(t_out, x_out) Cell_ode(t_out, x_out, param), tspan, x_initial, options);
[t_out, x_out, te, xe, ie] = ode15s(@(t_out, x_out) Cell_ode(t_out, x_out, param), tspan, x_initial, options);

r_out = transpose(x_initial);
u_out = r_out(2:end-1);
for t = 1:(param.t_duration / param.dt)
    % cathode
    [r_pi, u_pi] = finite_volume_update(r_out(t,1:param.Nr-1), u_out(t,1:param.Nr-2), param, ...
        param.delta_xp, param.Dsp_ref, param.F, param.A, param.Lp, param.a_sp);
    
    % anode
    [r_ni, u_ni] = finite_volume_update(r_out(t,param.Nr:end), u_out(t,param.Nr-1:end), param, ...
        param.delta_xn, param.Dsn_ref, param.F, param.A, param.Ln, param.a_sn);
    
    % combine 
    r_i = [r_pi, r_ni];
    u_i = [u_pi, u_ni];

    r_out = [r_out; r_i];
    u_out = [u_out; u_i];
    
end

%Transpose state matrix into row form to match established data structure 
x_out = r_out'; %x_out'; TODO!
I_dummy = param.I_data(1:size(x_out,2)); %Store exact current profiles used in I_dummy

if param.cycles == 0
%For no additional cycles, output the relevant variables
    param.I_data = I_dummy;
%     param.t_data = t_out; TODO!
    
else
    
%For any # of additional cycles > 0, continue cycling w/ alternating
%charge/discharge profiles

    for i = 1:param.cycles
        x_in2 = x_out(:,end);
        param.I_data = - param.I_data;
        event_formatted = @(t,x) physical_event_function(t,x,param); 
        options=odeset('RelTol',reltol,'AbsTol',abstol,'Events',event_formatted);
        [t_out2, x_out2,te2,xe2,ie2] = ode15s(@(t_out, x_out) Cell_ode(t_out, x_out, param), tspan, x_in2, options);

        %Transpose state matrix into row form to match established data structure 
        x_out2 = x_out2';
        %Shift t_out2 to account for duration of previous cycles
        t_out2 = t_out2+t_out(end)+1;

        x_out = [x_out x_out2];
        t_out = [t_out; t_out2];
        I_dummy = [I_dummy; param.I_data(1:size(x_out2,2))];
    end

%Pass concatenated time vector and current profiles to workspace in 'param'
    param.I_data = I_dummy;
    param.t_data = t_out;

end
%% Separate electrochemical state variables from x_out matrix
%Define solid concentrations
cs = x_out(1:(param.Nr-1)*param.Nc*2,:);            %All solid concentrations
cs_n = cs(1:(param.Nr-1)*param.Nc,:);               %Anode Concentrations
cs_p = cs((param.Nr-1)*param.Nc+1:end,:);           %Cathode Concentrations

for j = 1:length(cs)
    
    % Surface concentration -> surface stoichiometry
    for i = 1:param.Nc
        theta_surf_n(i,j) = cs_n(i*(param.Nr-1),j)/param.c_n_max;
        theta_surf_p(i,j) = cs_p(i*(param.Nr-1),j)/param.c_p_max;
    end
    
%% Open circuit potential and overpotential
 
    ocp_p(:,j) = U_p(theta_surf_p(:,j));
    ocp_n(:,j) = U_n(theta_surf_n(:,j));
    eta_p(:,j) = eta_cathode(theta_surf_p(:,j), param.ce0, param.T_cell, param.I_data(j), param, param.kp_ref);
    eta_n(:,j) = eta_anode(theta_surf_n(:,j), param.ce0, param.T_cell, param.I_data(j), param, param.kn_ref);
    

  % Cell voltages
  
    V_cell(:,j) = V_calculation(ocp_p(:,j),ocp_n(:,j),eta_p(:,j),eta_n(:,j),...
        param.I_data(j),param);
    V_oc = ocp_p - ocp_n;
    
end


%% Calculate SOC
[soc_bulk_n,soc_bulk_p] = soc_calculation(cs_n,cs_p,param);

end