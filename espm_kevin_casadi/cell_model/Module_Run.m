function [V_cell, R_l, T_core, T_surf,soc_bulk_n, soc_bulk_p, cs_n, cs_p,...
          ce_all,L_sei,i_s,Q,V_oc,R_el,R_sei,Csolv,aina_n,aina_p,c_sei,c_li,i_lpl,L_film,af_n,af_p,param, ...
          te, xe, ie] = Module_Run(x_initial,param)
%% Generate matrices for solid phase, electrolyte phase, thermal, and aging

import casadi.*


% Solid phase discretization matrices
[param.A_sd, param.B_sd_n, param.B_sd_p] = matrices_solidphase(param);

% Thermal model matrices (two-state thermal model)
[param.A_t, param.B_t] = matrices_thermal(param);

% Aging model matrices (one-state SEI layer growth model)
[param.A_sei, param.B_sei] = matrices_sei(param);

% Aging model matrices (one-state film layer growth model)
[param.A_film, param.B_film] = matrices_film(param);

% Aging model matrices (one-state LAM model)
[param.A_lam_n, param.A_lam_p, param.B_lam_n, param.B_lam_p] = matrices_lam(param);

% Electrolyte phase discretization matrices
[param.A_en_d, param.A_es_d, param.A_ep_d, param.B_e_d] = matrices_elecphase(param);

%% Separate compiled coefficient matrices (i.e. Define A1_en, A2_en, A3_en, etc.) 
run separate_elec_matrices.m

%% Solve ODEs 
tspan = param.t_data;
param.TIMEOLD = datetime('now','Format','HHmmss');




totN=size(x_initial,1);
sundials_u = SX.sym('sundials_u',totN);
input_crt=SX.sym('input_crt');

%================Define ODE using SUNDIALS

%================Define ODE using SUNDIALS
fun_ode=Module_ode_sundials(sundials_u,input_crt, param);
% dae = struct('x',sundials_u,'ode',[fun_ode],'alg',[]);
dae = struct('x',sundials_u,'p',input_crt,'ode',[fun_ode],'alg',[]);


if var(param.I_data)>1e-9



opts= struct('tf',1,'abstol',5.0e-8*0.001,'reltol',5.0e-8);  %设定仿真结束时间
F1 = integrator('F', 'idas', dae,opts);
u0=x_initial';


cc=1;
for kk=0:1:param.t_data(end)
%     kk
iapp=param.I_data(kk+1);
sol = F1('x0',u0,'p',iapp);
u0=full(sol.xf);

rex(:,cc)=full(sol.xf);
cc=cc+1;


%==========监测
if cc>10
% monitorxnow=rex(:,cc-1);
% monitorxpre=rex(:,cc-2);
monitorxnow=rex(:,cc-2);
monitorxpre=rex(:,cc-3);
value2now = max(monitorxnow(1:param.Nc*(param.Nr-1))) - (param.theta100_n*param.c_n_max-10);
value2pre = max(monitorxpre(1:param.Nc*(param.Nr-1))) - (param.theta100_n*param.c_n_max-10);

value3now = min(monitorxnow(param.Nc*(param.Nr-1)+1:2*param.Nc*(param.Nr-1))) - (param.theta100_p*param.c_p_max+10);
value3pre = min(monitorxpre(param.Nc*(param.Nr-1)+1:2*param.Nc*(param.Nr-1))) - (param.theta100_p*param.c_p_max+10);

if(value2now*value2pre<0)
    trigged_index=1;
    break
elseif(value3now*value3pre<0)
    trigged_index=2;
    break
end
end    
    
end


elseif var(param.I_data)<1e-9
    
tgrid=tspan;
opts = struct('tf',tspan(end),'grid',tgrid,'output_t0',true,'abstol',1e-6,'reltol',1e-6,'newton_scheme','gmres');  %设定仿真结束时间

% opts = struct('tf',tspan(end),'grid',tgrid,'output_t0',true,'abstol',1e-1,'reltol',1e-4,'newton_scheme','gmres');  %设定仿真结束时间
F1 = integrator('F', 'idas', dae,opts);
u0=x_initial';

sol = F1('x0',u0,'p',param.I_data(tspan(floor(tspan(end)/2))));


rex=full(sol.xf);



%====================检查浓度

% value2 = (max(rex(1:param.Nc*(param.Nr-1),:))) - (param.theta100_n*param.c_n_max-10);
% value3 = (min(rex(param.Nc*(param.Nr-1)+1:2*param.Nc*(param.Nr-1),:))) - (param.theta100_p*param.c_p_max+10);

monitorxnow=rex(:,1:end-1);
monitorxpre=rex(:,2:end);
value2now = max(monitorxnow(1:param.Nc*(param.Nr-1),:)) - (param.theta100_n*param.c_n_max-10);
value2pre = max(monitorxpre(1:param.Nc*(param.Nr-1),:)) - (param.theta100_n*param.c_n_max-10);

value3now = min(monitorxnow(param.Nc*(param.Nr-1)+1:2*param.Nc*(param.Nr-1),:)) - (param.theta100_p*param.c_p_max+10);
value3pre = min(monitorxpre(param.Nc*(param.Nr-1)+1:2*param.Nc*(param.Nr-1),:)) - (param.theta100_p*param.c_p_max+10);


index1=min(find(value2now.*value2pre<0));
index2=min(find(value3now.*value3pre<0));

% if (sum(index1)~=0 && sum(index2))~=0
%     useindex=min(index1,index2);
%     rex(:,useindex:end)=[];
% elseif (sum(index1)~=0 && sum(index2))==0
%     useindex=index1;
%     rex(:,useindex:end)=[];
% elseif (sum(index1)==0 && sum(index2))~=0
%     useindex=index2;
%     rex(:,useindex:end)=[];
% end


if (sum(index1)~=0 && sum(index2))~=0
    useindex=min(index1,index2);
    rex(:,useindex+2:end)=[];
elseif (sum(index1)~=0 && sum(index2))==0
    useindex=index1;
    rex(:,useindex+2:end)=[];
elseif (sum(index1)==0 && sum(index2))~=0
    useindex=index2;
    rex(:,useindex+2:end)=[];
end

end






x_out=rex;

num_for_cs=(param.Nr-1)*param.Nc*2;
num_for_ce=param.Nc*param.ce_states;
num_for_rest=2+7+param.Nc*param.Nsei;
tot_state_num=num_for_cs+num_for_ce+num_for_rest;


x_out=[x_out;300*ones(tot_state_num-(num_for_cs+num_for_ce),size(x_out,2))];


%============
if var(param.I_data)>1e-9
x_out=x_out(:,1:end-1);
x_add=[x_initial;300*ones(tot_state_num-(num_for_cs+num_for_ce),1)];
x_out=[x_add,x_out];
end

%=============

% %NB: Set to 5e-12 for HPPC; 5e-5 for capacity and UDDS
% reltol=5.0e-8; abstol=reltol*0.001;
% 
% % event_formatted = @(t,x) physical_event_function_tk(t,x,param); 
% % options=odeset('RelTol',reltol,'AbsTol',abstol,'Events',event_formatted);
% options=odeset('RelTol',reltol,'AbsTol',abstol);
% % [t_out, x_out,te,xe,ie] = ode15s(@(t_out, x_out) Module_ode(t_out, x_out, param), tspan, x_initial, options);
% [t_out, x_out] = ode15s(@(t_out, x_out) Module_ode(t_out, x_out, param), tspan, x_initial, options);
% 
% % options=odeset('RelTol',reltol,'AbsTol',abstol);
% % [t_out, x_out] = ode15s(@(t_out, x_out) Module_ode(t_out, x_out, param), tspan, x_initial, options);
% 
% % Transpose state matrix into row form to match established data structure 
% x_out = x_out';
% x_out=[x_out;zeros(17,size(x_out,2))];
I_dummy = param.I_data(1:size(x_out,2)); %Store exact current profiles used in I_dummy
t_out=0:1:param.t_data(end);
if param.cycles == 0
    % For no additional cycles, output the relevant variables
    param.I_data = I_dummy;
    param.t_data = t_out;
    
else
    % For any # of additional cycles > 0, continue cycling w/ alternating
    % charge/discharge profiles
    for i = 1:param.cycles
        x_in2 = x_out(:,end);
        param.I_data = - param.I_data;
        event_formatted = @(t,x) physical_event_function_tk(t,x,param); 
        options=odeset('RelTol',reltol,'AbsTol',abstol,'Events',event_formatted);
        [t_out2, x_out2,te2,xe2,ie2] = ode15s(@(t_out, x_out) Module_ode(t_out, x_out, param), tspan, x_in2, options);
%         [t_out2, x_out2,te2,xe2] = ode15s(@(t_out, x_out) Module_ode(t_out, x_out, param), tspan, x_in2, options);
        % Transpose state matrix into row form to match established data structure 
        x_out2 = x_out2';
        % Shift t_out2 to account for duration of previous cycles
        t_out2 = t_out2+t_out(end)+1;

        x_out = [x_out x_out2];
        t_out = [t_out; t_out2];
        I_dummy = [I_dummy; param.I_data(1:size(x_out2,2))];
    end

    % Pass concatenated time vector and current profiles to workspace in 'param'
    param.I_data = I_dummy;
    param.t_data = t_out;
end
%% Separate electrochemical, thermal & aging state variables from x_out matrix
%Define solid concentrations
cs = x_out(1:(param.Nr-1)*param.Nc*2,:);            %All solid concentrations
cs_n = cs(1:(param.Nr-1)*param.Nc,:);               %Anode Concentrations
cs_p = cs((param.Nr-1)*param.Nc+1:end,:);           %Cathode Concentrations
if ~isreal(cs_p)                                                    
    cs_p = abs(cs_p);
end

index_cs = (param.Nr-1)*param.Nc*2;                 %Index for final solid concentration
index_ce = index_cs + param.Nc*param.ce_states;     %Index for final electrolyte concentration

%Define electrolyte concentrations
ce = x_out(index_cs+1:index_ce,:);                  %All electrolyte concentrations
ce_n = ce(1:param.Nx_n*param.Nc,:);                 %Negative Electrolyte Region
ce_s = ce(param.Nx_n*param.Nc+1:(param.Nx_n + param.Nx_s)*param.Nc,:); %Separator Region
ce_p = ce((param.Nx_n + param.Nx_s)*param.Nc+1:param.ce_states*param.Nc,:); %Positive Region

%Define cell temperatures
T_cell = x_out(index_ce+1:index_ce+param.Nc*2,:);
index_thermal = index_ce+param.Nc*2;                %Index for final temperature state

%Define aging states
index_aging = index_thermal+7*param.Nc;                                %Index for final aging state
L_sei = x_out(index_thermal+1:index_thermal+param.Nc,:);               %SEI Layer Thickness
Q = x_out(index_thermal+param.Nc+1:index_thermal+2*param.Nc,:);        %ODE-updated Capacity 
aina_n = x_out(index_thermal+2*param.Nc+1:index_thermal+3*param.Nc,:); %Inactive area evolution, anode
aina_p = x_out(index_thermal+3*param.Nc+1:index_thermal+4*param.Nc,:); %Inactive area evolution, cahode
c_sei = x_out(index_thermal+4*param.Nc+1:index_thermal+5*param.Nc,:);  %SEI concentration
c_li = x_out(index_thermal+5*param.Nc+1:index_thermal+6*param.Nc,:);   %Plated li concentration
L_film = x_out(index_thermal+6*param.Nc+1:index_thermal+7*param.Nc,:); %Film Layer Thickness

%Define SEI Solvent Concentration States
Csolv = x_out(index_aging+1:end,:); %Solvent Concentration in SEI 



%% Temperature dependent transport and kinetics for each cell
for j = 1:length(cs)
    % Surface concentration -> surface stoichiometry
    for i = 1:param.Nc
        theta_surf_n(i,j) = cs_n(i*(param.Nr-1),j)/param.c_n_max;
        theta_surf_p(i,j) = cs_p(i*(param.Nr-1),j)/param.c_p_max;
        T_core(i,j) = T_cell(2*i-1,j);
        T_surf(i,j) = T_cell(2*i,j);
    end
    
    % Fracture area
    af_n(:,j) = param.a_sn*param.k_lam_n*(param.af_tinit+t_out(j))*param.flag_aina_n;
    af_p(:,j) = param.a_sp*param.k_lam_p*(param.af_tinit+t_out(j))*param.flag_aina_p;
    
    % Temperature dependence
    [Dsn(:,j), Dsp(:,j), kn(:,j), kp(:,j)] = arrhenius_temp(param,T_core(:,j));
    
    % SEI layer dependence
    [eps_el_n(:,j),R_sei(:,j)] = sei_lam_effect(L_film(:,j),param,aina_n(:,j),af_n(:,j));
    
    % LAM dependence
    eps_el_p(:,j) = lam_effect(param,aina_p(:,j),af_p(:,j));
    
    % Define electrolyte boundary concentrations for each cell
    for i = 1:param.Nc
        index1n = (i-1)*param.Nx_n+1;
        index2n = i*param.Nx_n;
        index1s = (i-1)*param.Nx_s+1;
        index2s = i*param.Nx_s;
        index1p = (i-1)*param.Nx_p+1;
        index2p = i*param.Nx_p;
      
        %Average Electrolyte Concentrations for conc. dependent parameters
        ce_n_avg(i,j) = mean(ce_n(index1n:index2n,j));
        ce_s_avg(i,j) = mean(ce_s(index1s:index2s,j));
        ce_p_avg(i,j) = mean(ce_p(index1p:index2p,j));
%         ce_all(:,j,i) = [ce_n(index1n:index2n,j); ce_s(index1s:index2s,j); ce_p(index1p:index2p,j)];
        
        %ce_all_avg(Row = 'Cell #', Column = 'Time point')
        ce_all_avg(i,j) = mean([ce_n(index1n:index2n,j); ce_s(index1s:index2s,j); ce_p(index1p:index2p,j)]);
        
    end
    
    % Electrolyte Conductivity
    [K_el_eff_n(:,j),K_el_eff_s(:,j),K_el_eff_p(:,j)] = conductivity_update(ce_n_avg(:,j)...
        ,ce_s_avg(:,j),ce_p_avg(:,j),eps_el_n(:,j),eps_el_p(:,j),T_core(:,j),param);
    
    % Open circuit potential and overpotential
    ocp_p(:,j) = U_p(theta_surf_p(:,j));
    ocp_n(:,j) = U_n(theta_surf_n(:,j));
    eta_p(:,j) = eta_cathode(theta_surf_p(:,j), ce_p_avg(:,j), T_core(:,j), param.I_data(j), param, kp(:,j), aina_p(:,j), af_p(:,j));
    eta_n(:,j) = eta_anode(theta_surf_n(:,j), ce_n_avg(:,j), T_core(:,j), param.I_data(j), param, kn(:,j), aina_n(:,j), af_n(:,j));  

    % Cell voltages
    [phi_e(:,j),R_el(:,j)] = electrolyte_potential(ce_all_avg(:,j), [ce_n(index1n:index2n,j); ce_s(index1s:index2s,j); ce_p(index1p:index2p,j)],K_el_eff_n(:,j),...
        K_el_eff_s(:,j), K_el_eff_p(:,j),T_core(:,j), param.I_data(j), param);
    [V_cell(:,j), R_l(:,j)] = V_calculation(ocp_p(:,j),ocp_n(:,j),eta_p(:,j),eta_n(:,j),...
        phi_e(:,j),R_sei(:,j),param.I_data(j),param.SOC_cc(j),param);
    V_oc = ocp_p - ocp_n;
    
    % Side reaction current density
    phi_sn(:,j) = ocp_n(:,j) + eta_n(:,j) + R_sei(:,j)*param.I_data(j);
    i_s(:,j) = side_current_solvent_aftersolve(R_sei(:,j), T_core(:,j), ...
        phi_sn(:,j), param, param.I_data(j), theta_surf_n(:,j),Csolv(:,j));
    i_lpl(:,j) = side_current_plating_aftersolve(R_sei(:,j),T_core(:,j),phi_sn(:,j),phi_e(:,j),param,param.I_data(j));
end
ce_all=[ce_n(index1n:index2n,:); ce_s(index1s:index2s,:); ce_p(index1p:index2p,:)];
% figure; hold on
% plot(eta_p); plot(-eta_n)

%% Calculate SOC
[soc_bulk_n,soc_bulk_p] = soc_calculation(cs_n,cs_p,param);

end