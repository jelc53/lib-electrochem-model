%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% SPM finite volume version %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Le Xu(lexu1209@stanford.edu) and Julian Edwin Lovett Cooper(jelc@stanford.edu)%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Date: 2022/12/01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% In case of help, feel free to contact the author %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%% addpath model
addpath cell_model

%% Load identification results

load([pwd '/opt_res/hist/x_opt_tk_v0_UDDS.mat'])


%% Load Experimental Data

% %=======================Condition 1: UDDS
% load('data/UDDS_cycle_1.mat')
% SOC_IC=0.8;
% use_time=3000;

%=======================Condition 2: Constant charge
use_time=3000;
randn('seed', 10);
I_data=-4.87/1*ones(1,use_time)+0*randn(1,use_time);
I_data=[I_data(1),I_data];
SOC_IC=0;

timeinterval=1;
t_data=0:timeinterval:use_time;
dt=1;
Q_IC=5;
Lsei_IC=0;
T_amb=23;
I_data=I_data(t_data+1);
V_data=ones(1,size(t_data,2));
SOC_cc=ones(1,size(t_data,2));

% % =======================Condition 4: Multiple charge and discharge cycles
% randn('seed',10);
% num_cycles=10;
% discharge_time1=500;
% discharge_time2=500;
% % use_time=discharge_time*2*num_cycles;
% use_time=(discharge_time1+discharge_time2)*num_cycles;
% I_data1=4.87/1*ones(1,discharge_time1); 
% I_data2=-4.87/1*ones(1,discharge_time2);
% I_data=repmat([I_data1,I_data2],1,num_cycles)+0*randn(1,use_time);
% I_data=[I_data(1),I_data];
% SOC_IC=1;  % if 1: discharge first then charge (also flip I_data signs)
% 
% timeinterval=1;
% t_data=0:timeinterval:use_time;
% dt=1;
% Q_IC=5;
% Lsei_IC=0;
% T_amb=23;
% I_data=I_data(t_data+1);
% V_data=ones(1,size(t_data,2));
% SOC_cc=ones(1,size(t_data,2));

%% Sim
%==================the last input indicates we are using FDM or FVM method
%===================0 represents FDM method
%===================1 represents FVM method

%===========================================================================
%=====================================FVM version
method_flag=1;
tic
[V_cell, R_l, T_core, T_surf,soc_bulk_n, soc_bulk_p, cs_n, cs_p,...
          ce_all,V_oc,R_el,R_sei,param] = ESPM_sim(x_opt,dt,t_data,I_data,SOC_cc,SOC_IC,Q_IC,Lsei_IC,T_amb,method_flag);     
caltime=toc

%====updata t_data
fvm_t_data=param.t_data;
%==================Date post processing
%=======================1 Calculate volume-average concentraion
%============The third input indicate cathode (0) or anode(1)
csp_ave=volume_ave_con_fvm(param,cs_p,0);
csn_ave=volume_ave_con_fvm(param,cs_n,1);

%=======================2 Calculate surface concentraion
%============The third input indicate cathode (0) or anode(1)
csp_surf=surface_con(param,cs_p,0);
csn_surf=surface_con(param,cs_n,1);


% ====================plot results
figure;
plot(fvm_t_data,V_cell,'Linewidth',2)
xlabel('Time (s)');
ylabel('Voltage (V)');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');


figure;
subplot(1,2,1)
plot(fvm_t_data,csp_ave,'r','Linewidth',2)
xlabel('Time (s)');
ylabel('Concentraion (mol/m^3)');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');

subplot(1,2,2)
plot(fvm_t_data,csn_ave,'r','Linewidth',2)
xlabel('Time (s)');
ylabel('Concentraion (mol/m^3)');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');


% Write results to file
all_data.t=fvm_t_data;
all_data.csp=cs_p;
all_data.csn=cs_n;
all_data.node_Nr=param.Nr;
all_data.node_Nxp=param.Nx_p;
all_data.node_Nxs=param.Nx_s;
all_data.node_Nxn=param.Nx_n;
all_data.V=V_cell;
all_data.cspave=csp_ave;
all_data.csnave=csn_ave;
all_data.cssurp=csp_surf;
all_data.cssurn=csn_surf;
save('matlab_results_fvm.mat','all_data') 

%============================================================================
%==================================================FDM version
method_flag=0;

tic
[V_cell, R_l, T_core, T_surf,soc_bulk_n, soc_bulk_p, cs_n, cs_p,...
          ce_all,V_oc,R_el,R_sei,param] = ESPM_sim(x_opt,dt,t_data,I_data,SOC_cc,SOC_IC,Q_IC,Lsei_IC,T_amb,method_flag);     
caltime=toc

%====updata t_data
fdm_t_data=param.t_data;
%==================Date post processing
%=======================1 Calculate volume-average concentraion
%============The third input indicate cathode (0) or anode(1)
[csp_ave,csn_ave]=volume_ave_con_fdm(param,cs_p,cs_n);


%=======================2 Calculate surface concentraion
%============The third input indicate cathode (0) or anode(1)
csp_surf=cs_p(end,:);
csn_surf=cs_n(end,:);


% ====================plot results
figure;
plot(V_cell,'Linewidth',2)
xlabel('Time (s)');
ylabel('Voltage (V)');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');


figure;
subplot(1,2,1)
plot(fdm_t_data,csp_ave,'r','Linewidth',2)
xlabel('Time (s)');
ylabel('Concentraion (mol/m^3)');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');

subplot(1,2,2)
plot(fdm_t_data,csn_ave,'r','Linewidth',2)
xlabel('Time (s)');
ylabel('Concentraion (mol/m^3)');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');


% Write results to file
all_data.t=fdm_t_data;
all_data.csp=cs_p;
all_data.csn=cs_n;
all_data.node_Nr=param.Nr;
all_data.node_Nxp=param.Nx_p;
all_data.node_Nxs=param.Nx_s;
all_data.node_Nxn=param.Nx_n;
all_data.V=V_cell;
all_data.cspave=csp_ave;
all_data.csnave=csn_ave;
all_data.cssurp=csp_surf;
all_data.cssurn=csn_surf;
save('matlab_results_fdm.mat','all_data') 


%% ===============Output model parameters and they will be used by Pybamm
SOC_cells=SOC_IC;
var.Rp=param.Rs_p;
var.Rn=param.Rs_n;

var.maximum_csp=param.c_p_max;
var.maximum_csn=param.c_n_max;
csn0 = (SOC_cells*(param.theta100_n - param.theta0_n) + param.theta0_n).*param.c_n_max ;
csp0 = (param.theta0_p - SOC_cells*(param.theta0_p - param.theta100_p)).*param.c_p_max; 
var.initial_csp=csp0;
var.initial_csn=csn0;

var.lp=param.Lp;
var.ls=param.Ls;
var.ln=param.Ln;


var.porosityp=param.eps_el_p;
var.porositys=param.eps_el_s;
var.porosityn=param.eps_el_n_ref;

var.epsp=param.epsilon_p;
var.epsn=param.epsilon_n;

var.Dsp=param.Dsp_ref;
var.Dsn=param.Dsn_ref;


var.kp=param.kp_ref;
var.kn=param.kn_ref;
var.A=param.A;
usevar=[888,var.Rp,var.Rn,var.maximum_csp,var.maximum_csn,var.initial_csp,var.initial_csn,...
    var.lp,var.ls,var.ln,var.porosityp,var.porositys,var.porosityn,var.epsp,var.epsn,var.Dsp,...
    var.Dsn,var.kp,var.kn,var.A]';
csvwrite('use_var.csv',usevar)



