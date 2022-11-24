%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ESPM validation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Gabriele Pozzato (gpozzato@stanford.edu) %%%%%% Date: 2020/12/29 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% In case of help, feel free to contact the author %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
global A_csn A_csp
global err_Up err_theta_surf_p
global trigged_index trigged_time
%% addpath model
addpath cell_model

%% Load identification results

load([pwd '/opt_res/hist/x_opt_tk_v0_UDDS.mat'])
% load([pwd '/opt_res/hist/x_opt_tk_v0_UDDS_1000_1s.mat'])

% % Set err_Up to 0 (remove effect of U_p identification)
% err_Up = 0;

%% Load Experimental Data
% Load cycling data for N<=25 preloaded for cell W8

% randn('seed',10);
% I_data=2.1*1*ones(1,size(I_data,2))+0.0001*randn(1,size(I_data,2));
% I_data=2.1*1*ones(1,size(I_data,2));
% SOC_IC=1;

% %=======================Condition 1: UDDS
% load('data/UDDS_cycle_1.mat')
% SOC_IC=0.8;
% use_time=3000;

%=======================Condition 2: Constant charge
use_time=3600;
randn('seed', 10);
I_data=-2.4/1*ones(1,use_time)+0.0001*randn(1,use_time);
SOC_IC=0;

% %=======================Condition 3: Constant charge 1/20 C
% use_time=3600*20;
% randn('seed',10);
% I_data=-2.4/20*ones(1,use_time)+0.0001*randn(1,use_time);
% SOC_IC=0;

% =======================Condition 4: Multiple charge and discharge cycles
% randn('seed',10);
% num_cycles=10;
% discharge_time=3600;
% use_time=discharge_time*2*num_cycles;
% I_data=-2.4/1*ones(1,discharge_time);
% I_data=repmat([I_data,-I_data],1,num_cycles)+0.0001*randn(1,use_time);
% SOC_IC=0;

timeinterval=1;
t_data=0:timeinterval:use_time-1;
dt=1;
Q_IC=5;
Lsei_IC=0;
T_amb=23;
% I_data=I_data(1:use_time);
% V_data=ones(1,use_time);
% SOC_cc=ones(1,use_time);
I_data=I_data(t_data+1);
V_data=ones(1,size(t_data,2));
SOC_cc=ones(1,size(t_data,2));

%% Sim
% NB
% Set event function max time to 1000s
% Set relative tolerance to 5e-8

% T_amb=30
tic
% [V_cell, R_l, T_core, T_surf,soc_bulk_n, soc_bulk_p, cs_n, cs_p,...
%  ce_all,L_sei,i_s,Q,V_oc,R_el,R_sei,Csolv,aina_n,aina_p,...
%  c_sei,c_li,i_lpl,L_film,af_n,af_p,param] = ESPM_sim(x_opt,dt,t_data,I_data,SOC_cc,SOC_IC,Q_IC,Lsei_IC,T_amb);
[V_cell, R_l, T_core, T_surf,soc_bulk_n, soc_bulk_p, cs_n, cs_p,...
          ce_all,V_oc,R_el,R_sei,param] = ESPM_sim(x_opt,dt,t_data,I_data,SOC_cc,SOC_IC,Q_IC,Lsei_IC,T_amb);
[w csp_avg] = radial_average(cs_p,param.Nr,param.Rs_p,param.delta_xp);
[w csn_avg] = radial_average(cs_n,param.Nr,param.Rs_n,param.delta_xn);
caltime=toc

%% Sense check plots
plot(V_cell)

figure;
plot(cs_n(end,:))

figure;
plot(cs_p(end,:))

% figure;
% plot(csp_avg)
% 
% figure;
% plot(csn_avg)

%% Write to file
all_data.ce=ce_all;
all_data.csp=cs_p;
all_data.csn=cs_n;
all_data.node_Nr=param.Nr;
all_data.node_Nxp=param.Nx_p;
all_data.node_Nxs=param.Nx_s;
all_data.node_Nxn=param.Nx_n;
all_data.par=param;
all_data.V=V_cell;
all_data.caltime=caltime;

save('Oldsolver_re.mat','all_data') %'Oldsolver_re.mat'
