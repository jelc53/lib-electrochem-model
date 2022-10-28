function [dx_dt] = Module_ode(t_in, x_in, param)
global A_csn A_csp
% Solving solid phase PDE using ODE solver
% Interpolate current


num_for_cs=(param.Nr-1)*param.Nc*2;
num_for_ce=param.Nc*param.ce_states;
num_for_rest=2+7+param.Nc*param.Nsei;
tot_state_num=num_for_cs+num_for_ce+num_for_rest;

x_in=[x_in;300*ones(tot_state_num-(num_for_cs+num_for_ce),1)];

input_crt = interp1(param.t_data, param.I_data, t_in);   

% Interpolate Coulomb-counted SOC
input_SOC_cc = interp1(param.t_data, param.SOC_cc, t_in);

%% New: update diffusivity (v2 - Release 2)
% Imported from March 2020 - Module Code: Solv_diff_update.m
% Compute Dsolv at time t as a function of C-rate (input current)
param.Dsolv_ref = Solv_Diff_update(input_crt);

%% Separate electrochemical, thermal & aging state variables from x_in vector
%Define solid concentrations
cs = x_in(1:(param.Nr-1)*param.Nc*2);           %Solid Concentration states
cs_n = cs(1:(param.Nr-1)*param.Nc);             %Anode Solid Conc.
cs_p = cs((param.Nr-1)*param.Nc+1:end);         %Cathode Solid Conc. 
index_cs = (param.Nr-1)*param.Nc*2;             %Index for final solid concentration
if ~isreal(cs_p)                                %Until the event is triggered
    cs_p = abs(cs_p);
end

%Define electrolyte concentrations
index_ce = index_cs + param.Nc*param.ce_states; %Index for final electrolyte concentration
ce = x_in(index_cs+1:index_ce);                 %All electrolyte concentrations
ce_n = ce(1:param.Nx_n*param.Nc);               %Negative Electrolyte Region
ce_s = ce(param.Nx_n*param.Nc+1:(param.Nx_n + param.Nx_s)*param.Nc); %Separator Region
ce_p = ce((param.Nx_n + param.Nx_s)*param.Nc+1:param.ce_states*param.Nc); %Positive Region

%Define cell temperatures
T_cell = x_in(index_ce+1:index_ce+param.Nc*2);  %Cell temperatures (surface + core)
index_thermal = index_ce+param.Nc*2;            %Index for final temperature state

%Define aging states
index_aging = index_thermal+7*param.Nc;                                %Index for final aging state
L_sei = x_in(index_thermal+1:index_thermal+param.Nc);                  %SEI Layer Thickness
Q = x_in(index_thermal+param.Nc+1:index_thermal+2*param.Nc);           %ODE-updated Capacity   
aina_n = x_in(index_thermal+2*param.Nc+1:index_thermal+3*param.Nc,:);  %Inactive area evolution, anode
aina_p = x_in(index_thermal+3*param.Nc+1:index_thermal+4*param.Nc,:);  %Inactive area evolution, cahode
c_sei = x_in(index_thermal+4*param.Nc+1:index_thermal+5*param.Nc,:);   %SEI concentration
c_li = x_in(index_thermal+5*param.Nc+1:index_thermal+6*param.Nc,:);    %Plated li concentration
L_film = x_in(index_thermal+6*param.Nc+1:index_thermal+7*param.Nc,:);  %Film Layer Thickness

%Define SEI Solvent Concentration States
Csolv = x_in(index_aging+1:end); %Solvent Concentration in SEI 

%% Temperature dependent transport and kinetics for each cell
%Allocate memory
theta_surf_n(param.Nc,1) = 0;
theta_surf_p(param.Nc,1) = 0;
T_core(param.Nc,1) = 0;
T_surf(param.Nc,1) = 0;

for i = 1:param.Nc
    % Surface concentration -> surface stoichiometry
    theta_surf_n(i,:) = cs_n(i*(param.Nr-1))/param.c_n_max;
    theta_surf_p(i,:) = cs_p(i*(param.Nr-1))/param.c_p_max;  
    T_core(i,:) = T_cell(2*i-1);
    T_surf(i,:) = T_cell(2*i);
end

%% Update parameters based on aging & temperature
% LAM
af_n = param.a_sn*param.k_lam_n*(param.af_tinit+t_in)*param.flag_aina_n;
af_p = param.a_sp*param.k_lam_p*(param.af_tinit+t_in)*param.flag_aina_p;

% Temperature dependence
[Dsn, Dsp, kn, kp, Dsolv] = arrhenius_temp(param,T_core);

% SEI layer dependence
[eps_el_n, R_sei] = sei_lam_effect(L_film,param,aina_n,af_n);

% LAM dependence
eps_el_p = lam_effect(param,aina_p,af_p);

%% Calculate average electrolyte concentrations
%Allocate Memory
ce_n_avg(param.Nc,1) = 0;
ce_s_avg(param.Nc,1) = 0;
ce_p_avg(param.Nc,1) = 0;
ce_all((param.Nx_n+param.Nx_s+param.Nx_p),param.Nc) = 0;
ce_all_avg(param.Nc,1) = 0;

for i = 1:param.Nc
    index1n = (i-1)*param.Nx_n+1; %1st negative electrolyte grid point for cell 'i'
    index2n = i*param.Nx_n;       %Last negative electrolyte grid point for cell 'i'
    index1s = (i-1)*param.Nx_s+1; %1st separator grid point for cell 'i'
    index2s = i*param.Nx_s;       %Last separator grid point for cell 'i'
    index1p = (i-1)*param.Nx_p+1; %1st positive electrolyte grid point for cell 'i'
    index2p = i*param.Nx_p;       %Last positive electrolyte grid point for cell 'i'
    
   %Average Electrolyte Concentrations
    ce_n_avg(i,1) = mean(ce_n(index1n:index2n));
    ce_s_avg(i,1) = mean(ce_s(index1s:index2s));
    ce_p_avg(i,1) = mean(ce_p(index1p:index2p));
    ce_all(:,i) = [ce_n(index1n:index2n); ce_s(index1s:index2s); ce_p(index1p:index2p)];
    ce_all_avg(i,1) = mean(ce_all(:,i));
    
end

%% Update Electrolyte Concentration- / Temperature- / Aging - dependent parameters
% Electrolyte Diffusivity

[De_n,De_s,De_p] = diffusivity_update(ce_n,ce_s,ce_p,T_core,param);

    De_n_minus0_5 = zeros(param.Nc*param.Nx_n,1);
    De_n_plus0_5 = zeros(param.Nc*param.Nx_n,1);
    
    De_s_minus0_5 = zeros(param.Nc*param.Nx_s,1);
    De_s_plus0_5 = zeros(param.Nc*param.Nx_s,1);
    
    De_p_minus0_5 = zeros(param.Nc*param.Nx_p,1);
    De_p_plus0_5 = zeros(param.Nc*param.Nx_p,1);

for i = 1:param.Nc
    
    index1n = (i-1)*param.Nx_n+1;
    index2n = i*param.Nx_n;       
    index1s = (i-1)*param.Nx_s+1; 
    index2s = i*param.Nx_s;       
    index1p = (i-1)*param.Nx_p+1; 
    index2p = i*param.Nx_p;       
    
    %Create vectors of diffusivity at grid points i - 1 or i + 1
    De_n_minus1 = [0; De_n(index1n:index2n-1,1)];
    De_n_plus1 = [De_n(index1n+1:index2n,1); 0];
    De_n_dummy = De_n(index1n:index2n,1);
    
    %Calculate diffusivities at grid points i + 1/2 and i - 1/2
    De_n_minus0_5(index1n:index2n) = 2*De_n_dummy.*De_n_minus1./(De_n_minus1 + De_n_dummy);
    De_n_plus0_5(index1n:index2n) = 2*De_n_dummy.*De_n_plus1./(De_n_plus1 + De_n_dummy);
    
    De_s_minus1 = [0; De_s(index1s:index2s-1,1)];
    De_s_plus1 = [De_s(index1s+1:index2s,1); 0];
    De_s_dummy = De_s(index1s:index2s,1);
    
    De_s_minus0_5(index1s:index2s) = 2*De_s_dummy.*De_s_minus1./(De_s_minus1 + De_s_dummy);
    De_s_plus0_5(index1s:index2s) = 2*De_s_dummy.*De_s_plus1./(De_s_plus1 + De_s_dummy);
    
    De_p_minus1 = [0; De_p(index1p:index2p-1,1)];
    De_p_plus1 = [De_p(index1p+1:index2p,1); 0];
    De_p_dummy = De_p(index1p:index2p,1);
    
    De_p_minus0_5(index1p:index2p) = 2*De_p_dummy.*De_p_minus1./(De_p_minus1 + De_p_dummy);
    De_p_plus0_5(index1p:index2p) = 2*De_p_dummy.*De_p_plus1./(De_p_plus1 + De_p_dummy);
    
    %Define Electrolyte Interface Effective Diffusivities
    %Note: these formulas assume same geometry for each cell
    %Note: interface diffusivities are defined separately to simplify
    %indexing
    beta_n_s = param.delta_n/(param.delta_n + param.delta_s);
    De_eff_n_s(i) = (De_n(index2n)*eps_el_n(i)^param.brugg_n)*(De_s(index1s)*param.eps_el_s^param.brugg_s)/(beta_n_s*...
        (De_s(index1s)*param.eps_el_s^param.brugg_s) + (1-beta_n_s)*(De_n(index2n)*eps_el_n(i)^param.brugg_n));
    param.delta_n_s = 1/2*(param.delta_n + param.delta_s);
    
    beta_s_p = param.delta_s/(param.delta_s + param.delta_p);
    De_eff_s_p(i) = (De_s(index2s)*param.eps_el_s^param.brugg_s)*(De_p(index1p)*eps_el_p(i)^param.brugg_p)/(beta_s_p*...
        (De_p(index1p)*eps_el_p(i)^param.brugg_p) + (1-beta_s_p)*(De_s(index2s)*param.eps_el_s^param.brugg_s));
    param.delta_s_p = 1/2*(param.delta_s + param.delta_p);
    
end

% Electrolyte Conductivity
[K_el_eff_n,K_el_eff_s,K_el_eff_p] = conductivity_update(ce_n_avg,ce_s_avg,ce_p_avg,eps_el_n,eps_el_p,T_core,param);

% Open circuit potential and overpotential
ocp_p = U_p(theta_surf_p);
ocp_n = U_n(theta_surf_n);
eta_p = eta_cathode(theta_surf_p, ce_p_avg, T_core, input_crt, param, kp, aina_p, af_p);
eta_n = eta_anode(theta_surf_n, ce_n_avg, T_core, input_crt, param, kn, aina_n, af_n);

% Cell voltages
[phi_e,R_el] = electrolyte_potential(ce_all_avg, ce_all,K_el_eff_n,K_el_eff_s,K_el_eff_p,T_core, input_crt, param);

% TODO: Change this?? --> not needed since we are not using the thermal
% model... but may cause issues later once we need it
% therefore TODO: Check that this actually works
V_cell = ocp_p - ocp_n + eta_p - eta_n + phi_e - param.R_l*input_crt - R_sei*input_crt;
% [V_cell, R_l] = V_calculation(ocp_p,ocp_n,eta_p,eta_n,phi_e,R_sei,input_crt,input_SOC_cc, param);

% Heat generation in each cell
V_loss = (ocp_p - ocp_n - V_cell); 
B_t = param.B_t;
for i = 1:param.Nc
   B_t(2*i-1,1) = param.B_t(2*i-1,1)*V_loss(i,:); %Incorporate heat loss voltage into B matrix
end

% Side reaction current density
% phi_sn = ocp_n + eta_n + L_sei*input_crt/param.kappa_sei;
phi_sn = ocp_n + eta_n + R_sei*input_crt;
i_s = side_current_solvent(R_sei,T_core,phi_sn,param,input_crt,theta_surf_n,Csolv);
i_lpl = side_current_plating(R_sei,T_core,phi_sn,phi_e,param,input_crt);

%% Solve solid phase ODEs

% Coefficients of discretized ODEs
alpha_n = Dsn/(param.delta_xn^2);
alpha_p = Dsp/(param.delta_xp^2);

W_mat_n = zeros(param.Nc*(param.Nr-1)); %Coeff. matrix for all anode solid conc. states for Nc cells
W_mat_p = zeros(param.Nc*(param.Nr-1)); %Coeff. matrix for all cathode solid conc. states for Nc cells

for i = 1:param.Nc
    %Define W_mat for all solid conc. states for each cell
    %Repmat creates W_mat by tiling Alpha(i) in (Nr - 1)x(Nr -1) matrix
    W_mat_n((i-1)*(param.Nr - 1)+1:i*(param.Nr - 1),(i-1)*(param.Nr - 1)+1:i*(param.Nr - 1)) = ...
        repmat(alpha_n(i,:),(param.Nr - 1),(param.Nr - 1));
    W_mat_p((i-1)*(param.Nr - 1)+1:i*(param.Nr - 1),(i-1)*(param.Nr - 1)+1:i*(param.Nr - 1)) = ...
        repmat(alpha_p(i,:),(param.Nr - 1),(param.Nr - 1));    
end

A_csn = param.A_sd.*W_mat_n; %Full Anode Solid Conc. State A-matrix
A_csp = param.A_sd.*W_mat_p; %%Full Cathode Solid Conc. State A-matrix

% Create i_s vector with proper dimensions
i_s_dummy = ones(param.Nc*(param.Nr-1),1);
i_lpl_dummy = ones(param.Nc*(param.Nr-1),1);
for i = 1:param.Nc
    i_s_dummy((i-1)*(param.Nr-1)+1:i*(param.Nr-1),1) = i_s_dummy((i-1)*...
        (param.Nr-1)+1:i*(param.Nr-1),1)*i_s(i);
    i_lpl_dummy((i-1)*(param.Nr-1)+1:i*(param.Nr-1),1) = i_lpl_dummy((i-1)*...
        (param.Nr-1)+1:i*(param.Nr-1),1)*i_lpl(i);
end

%Use modified intercalation current for anode
%NOTE: Sign convention is such that '+ i_s' is correct sign in
%inter_crt expression
dcsn_dt = A_csn*cs_n + param.B_sd_n.*(input_crt/(param.a_sn+af_n-aina_n)-(i_s_dummy+i_lpl_dummy)*param.A*param.Ln);
dcsp_dt = A_csp*cs_p + param.B_sd_p*input_crt/(param.a_sp+af_p-aina_p);
dcs_dt = [dcsn_dt; dcsp_dt];

%% Solve Electrolyte Phase ODE
% Coefficients of discretized electrolyte ODE's
alpha_el_n = eps_el_n.^(param.brugg_n-1)./(param.delta_n^2);

alpha_el_s = param.eps_el_s^(param.brugg_s-1)./(param.delta_s^2);

alpha_el_p = eps_el_p.^(param.brugg_p-1)./(param.delta_p^2);

beta_el_n = (1-param.t0)./(eps_el_n*param.F*param.A*param.Ln);
beta_el_s = 0; 
beta_el_p = (1-param.t0)/(eps_el_p*param.F*param.A*param.Lp); 
%NOTE: +/- for anode / cathode is incorporated later

%DEFINE W_mat matrices
W_e_n = zeros(param.Nc*param.Nx_n);
W_e_s = zeros(param.Nc*param.Nx_s);
W_e_p = zeros(param.Nc*param.Nx_p);

for i = 1:param.Nc; 
    %Define W_mat for all electrolyte conc. states for each cell
    %Repmat creates W_mat by tiling Alpha(i) in (Nr - 1)x(Nr -1) matrix
    W_e_n((i-1)*param.Nx_n+1:i*param.Nx_n,(i-1)*param.Nx_n+1:i*param.Nx_n) = ...
        repmat(alpha_el_n(i,:),param.Nx_n,param.Nx_n);

    W_e_s((i-1)*param.Nx_s+1:i*param.Nx_s,(i-1)*param.Nx_s+1:i*param.Nx_s) = ...
        repmat(alpha_el_s,param.Nx_s,param.Nx_s);
    
    W_e_p((i-1)*param.Nx_p+1:i*param.Nx_p,(i-1)*param.Nx_p+1:i*param.Nx_p) = ...
        repmat(alpha_el_p,param.Nx_p,param.Nx_p);
    
end

 Y_e_n = [];
 Y_e_s = [];
 Y_e_p = [];
 
    for i = 1:param.Nc
    n_dummy = beta_el_n(i)*ones(param.Nx_n,1);
    Y_e_n = [Y_e_n; n_dummy];
    
    s_dummy = beta_el_s*ones(param.Nx_s,1);
    Y_e_s = [Y_e_s; s_dummy];
    
    p_dummy = beta_el_p*ones(param.Nx_p,1);
    Y_e_p = [Y_e_p; p_dummy];
    end

    A_tot_n = De_n_minus0_5.*param.A1_en_d + De_n_plus0_5.*param.A2_en_d;
    A_tot_s = De_s_minus0_5.*param.A1_es_d + De_s_plus0_5.*param.A2_es_d;
    A_tot_p = De_p_minus0_5.*param.A1_ep_d + De_p_plus0_5.*param.A2_ep_d;
    
A_e_n = A_tot_n.*W_e_n;

B_e_n = Y_e_n.*param.B_en_d;

A_e_s = A_tot_s.*W_e_s;
B_e_s = Y_e_s.*param.B_es_d;

A_e_p = A_tot_p.*W_e_p;

B_e_p = Y_e_p.*param.B_ep_d;
    
%% Define Electrolyte Interface ODE's
%NOTE: interfaces between electrolyte regions are handled separately to simplify indexing
%Can try streamlining to handle interfaces with same algorithm as bulk in future versions
for i = 1:param.Nc
    index2n = i*param.Nx_n;
    index1s = (i-1)*param.Nx_s+1;
    index2s = i*param.Nx_s;
    index1p = (i-1)*param.Nx_p+1;
    
dCen_interf_ns(i) = 1/eps_el_n(i)*De_eff_n_s(i)/(param.delta_n*param.delta_n_s)*(ce_s(index1s) - ce_n(index2n))...
    -De_n_minus0_5(index2n)*eps_el_n(i)^(param.brugg_n-1)/(param.delta_n^2)*(ce_n(index2n) - ce_n(index2n-1)) + beta_el_n(i)*input_crt;

dCes_interf_ns(i) = De_s_plus0_5(index1s)*param.eps_el_s^(param.brugg_s-1)/(param.delta_s^2)*(ce_s(index1s+1) - ce_s(index1s))...
    -1/param.eps_el_s*De_eff_n_s(i)/(param.delta_s*param.delta_n_s)*(ce_s(index1s) - ce_n(index2n));


dCes_interf_sp(i) = 1/param.eps_el_s*De_eff_s_p(i)/(param.delta_s*param.delta_s_p)*(ce_p(index1p) - ce_s(index2s))...
    -De_s_minus0_5(index2s)*param.eps_el_s^(param.brugg_s-1)/(param.delta_s^2)*(ce_s(index2s) - ce_s(index2s-1));

dCep_interf_sp(i) = De_p_plus0_5(index1p)*eps_el_p(i)^(param.brugg_s-1)/(param.delta_p^2)*(ce_p(index1p+1) - ce_p(index1p))...
    -1/eps_el_p(i)*De_eff_s_p(i)/(param.delta_p*param.delta_s_p)*(ce_p(index1p) - ce_s(index2s)) - beta_el_p*input_crt;
end

%%   
% Define Electrolyte phase ODE's
%NOTE: 
%Input current term is added in negative electrolyte region
%Input current term is subtracted in positive electrolyte region
dce_n = A_e_n*ce_n + B_e_n*input_crt;

dce_s = A_e_s*ce_s + B_e_s*input_crt;

dce_p = A_e_p*ce_p - B_e_p*input_crt;

for i = 1:param.Nc
    
    index2n = i*param.Nx_n;
    index1s = (i-1)*param.Nx_s+1;
    index2s = i*param.Nx_s;
    index1p = (i-1)*param.Nx_p+1;
    
    dce_n(index2n) = dCen_interf_ns(i);
    dce_s(index1s) = dCes_interf_ns(i);
    dce_s(index2s) = dCes_interf_sp(i);
    dce_p(index1p) = dCep_interf_sp(i);
end

dce_dt = [dce_n;dce_s;dce_p];

%% Solve 1st-order thermal ODE
% dTcell_dt = param.A_t*T_cell + B_t*[input_crt; param.T_amb];
dTcell_dt = [0; 0];

%% Solve aging ODE (SEI -- diffusion-limited; plating; loss of active material (LAM); capacity loss) [Pozzato CCTA 2021]

% TODO: set to zeros
% 
% % SEI and plating
% dLsei_dt = param.A_sei*L_sei + param.B_sei*[i_s i_lpl]';
% dLfilm_dt = param.A_film*L_film + param.B_film*[i_s i_lpl]';
% dcsei_dt = -(param.a_sn+af_n-aina_n)*(i_s+i_lpl*param.beta_li2sei)/(2*param.F);
% dcli_dt = -(param.a_sn+af_n-aina_n)*i_lpl*(1-param.beta_li2sei);
% 
% % LAM
% daina_dt_n = (param.A_lam_n*aina_n + param.B_lam_n*[1 param.af_tinit+t_in]')*param.flag_aina_n;
% daina_dt_p = (param.A_lam_p*aina_p + param.B_lam_p*[1 param.af_tinit+t_in]')*param.flag_aina_p;
% 
% % Q loss
% dQ_dt = (i_s+i_lpl)*((param.a_sn+af_n-aina_n)*param.Ln*param.A)/3600;


% SEI and plating
dLsei_dt = 0;
dLfilm_dt = 0;
dcsei_dt = 0;
dcli_dt = 0;

% LAM
daina_dt_n = 0;
daina_dt_p = 0;

% Q loss
dQ_dt = 0;

%% Solvent Diffusion Dynamics in SEI Layer [Prada 2013]
dCsolv_dt(param.Nsei*param.Nc,1) = 0; %Allocate memory

for j = 1:param.Nc
    for i = 1:param.Nsei
        index = (j-1)*param.Nsei;
        xi = param.delta_xi*(i-1);

        if i == 1
            dCsolv_dt(index+i,1) = (Dsolv(j)/((L_film(j)*param.delta_xi)^2)*(2*Csolv(index+i+1)-2*Csolv(index+i))+...
                (2/(L_film(j)*param.delta_xi) + 1/Dsolv(j)*dLfilm_dt(j))*((i_s(j)+i_lpl(j)*param.beta_li2sei)/param.F...
                - dLfilm_dt(j)*Csolv(index+i)));
        elseif i == param.Nsei
            dCsolv_dt(index+i,1) = 0;
        else
            dCsolv_dt(index+i,1) = (Dsolv(j)/((L_film(j)*param.delta_xi)^2)*(Csolv(index+i+1) - 2*Csolv(index+i)+...
                Csolv(index+i-1)) + (xi-1)*(1/L_film(j))*dLfilm_dt(j)*(Csolv(index+i+1)-Csolv(index+i-1))/(2*param.delta_xi));
        end
    end
end

%% Final formulation
% dx_dt = [dcs_dt; dce_dt; dTcell_dt; dLsei_dt; dQ_dt; daina_dt_n; daina_dt_p; dcsei_dt; dcli_dt; dLfilm_dt; dCsolv_dt];
dx_dt = [dcs_dt;dce_dt];
end