clear all
clc

N_DIMS_LIST = [5 10 20 40 60 80];

sz = [length(N_DIMS_LIST) 3];
varNames = {'pybamm', 'fvm', 'fdm'};
varTypes = ["double", "double", "double"];

V_rmse_table=table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
cspsur_rmse_table=table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
csnsur_rmse_table=table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
cspave_rmse_table=table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
csnave_rmse_table=table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

%%%%%%%%%%%%%% REFERENCE %%%%%%%%%%%%%%%%%%%%%
all_matlab_ref=load('matlab_ref_fdm_nr100.mat');
all_matlab_ref=all_matlab_ref.all_data;

t_ref=all_matlab_ref.t;
t_ref=t_ref';

V_ref=all_matlab_ref.V;
V_ref=V_ref';

cspave_ref=all_matlab_ref.cspave;
csnave_ref=all_matlab_ref.csnave;
cspave_ref=cspave_ref';
csnave_ref=csnave_ref';

cspsur_ref=all_matlab_ref.cssurp;
csnsur_ref=all_matlab_ref.cssurn;
cspsur_ref=cspsur_ref';
csnsur_ref=csnsur_ref';

%===check if the last timestep is float
if t_ref(end)-round(t_ref(end))~=0
%===the last timestep is float, delete it    
t_ref(end)=[];
V_ref(end)=[];   
cspave_ref(end)=[];  
csnave_ref(end)=[];  
cspsur_ref(end)=[];  
csnsur_ref(end)=[];  
end



%%%%%%%%%%%%%% FINITE DIFFERENCE %%%%%%%%%%%%%
[V_arr, cspsur_arr, csnsur_arr, cspave_arr, csnave_arr] = deal([], [], [], [], []);

for i = 1:length(N_DIMS_LIST)
    all_matlab_fdm=load(strcat('matlab_results_fdm_nr',string(N_DIMS_LIST(i)),'.mat'));
    all_matlab_fdm=all_matlab_fdm.all_data;
    
    t_fdm=all_matlab_fdm.t;
    t_fdm=t_fdm';
    
    V_fdm=all_matlab_fdm.V;
    V_fdm=V_fdm';
    
    cspave_fdm=all_matlab_fdm.cspave;
    csnave_fdm=all_matlab_fdm.csnave;
    cspave_fdm=cspave_fdm';
    csnave_fdm=csnave_fdm';
    
    cspsur_fdm=all_matlab_fdm.cssurp;
    csnsur_fdm=all_matlab_fdm.cssurn;
    cspsur_fdm=cspsur_fdm';
    csnsur_fdm=csnsur_fdm';
    
    %===check if the last timestep is float
    if t_fdm(end)-round(t_fdm(end))~=0
    %===the last timestep is float, delete it    
    t_fdm(end)=[];
    V_fdm(end)=[];   
    cspave_fdm(end)=[];  
    csnave_fdm(end)=[];  
    cspsur_fdm(end)=[];  
    csnsur_fdm(end)=[];  
    end

    V_arr=[V_arr; rms(V_fdm-V_ref)*1000];
    cspsur_arr=[cspsur_arr; rms(cspsur_fdm-cspsur_ref)];
    csnsur_arr=[csnsur_arr; rms(csnsur_fdm-csnsur_ref)];
    cspave_arr=[cspave_arr; rms(cspave_fdm-cspave_ref)];
    csnave_arr=[csnave_arr; rms(csnave_fdm-csnave_ref)];
end
V_rmse_table.fdm = V_arr;
cspsur_rmse_table.fdm = cspsur_arr;
csnsur_rmse_table.fdm = csnsur_arr;
cspave_rmse_table.fdm = cspave_arr;
csnave_rmse_table.fdm = csnave_arr;



%%%%%%%%%%%%%% FINITE VOLUME %%%%%%%%%%%%%%%%
[V_arr, cspsur_arr, csnsur_arr, cspave_arr, csnave_arr] = deal([], [], [], [], []);

for i = 1:length(N_DIMS_LIST)
    all_matlab_fvm=load(strcat('matlab_results_fvm_nr',string(N_DIMS_LIST(i)),'.mat'));
    all_matlab_fvm=all_matlab_fvm.all_data;
    
    t_fvm=all_matlab_fvm.t;
    t_fvm=t_fvm';
    
    V_fvm=all_matlab_fvm.V;
    V_fvm=V_fvm';
    
    cspave_fvm=all_matlab_fvm.cspave;
    csnave_fvm=all_matlab_fvm.csnave;
    cspave_fvm=cspave_fvm';
    csnave_fvm=csnave_fvm';
    
    cspsur_fvm=all_matlab_fvm.cssurp;
    csnsur_fvm=all_matlab_fvm.cssurn;
    cspsur_fvm=cspsur_fvm';
    csnsur_fvm=csnsur_fvm';
    
    %===check if the last timestep is float
    if t_fvm(end)-round(t_fvm(end))~=0
    %===the last timestep is float, delete it    
    t_fvm(end)=[];
    V_fvm(end)=[];   
    cspave_fvm(end)=[];  
    csnave_fvm(end)=[];  
    cspsur_fvm(end)=[];  
    csnsur_fvm(end)=[];  
    end

    V_arr=[V_arr; rms(V_fvm-V_ref)*1000];
    cspsur_arr=[cspsur_arr; rms(cspsur_fvm-cspsur_ref)];
    csnsur_arr=[csnsur_arr; rms(csnsur_fvm-csnsur_ref)];
    cspave_arr=[cspave_arr; rms(cspave_fvm-cspave_ref)];
    csnave_arr=[csnave_arr; rms(csnave_fvm-csnave_ref)];
end
V_rmse_table.fvm = V_arr;
cspsur_rmse_table.fvm = cspsur_arr;
csnsur_rmse_table.fvm = csnsur_arr;
cspave_rmse_table.fvm = cspave_arr;
csnave_rmse_table.fvm = csnave_arr;



%%%%%%%%%%%%%%%%% PYBAMM %%%%%%%%%%%%%%%%%%%%%
[V_arr, cspsur_arr, csnsur_arr, cspave_arr, csnave_arr] = deal([], [], [], [], []);

t_pybamm=load('pybamm_t_single.txt');
V_pybamm=load('pybamm_V_single.txt');

cspave_pybamm=load('pybamm_cspave_single.txt');
csnave_pybamm=load('pybamm_csnave_single.txt');

cspsur_pybamm=load('pybamm_cspsur_single.txt');
csnsur_pybamm=load('pybamm_csnsur_single.txt');

V_rmse_table.pybamm = rms(V_pybamm-V_ref)'*1000;
cspsur_rmse_table.pybamm = rms(cspsur_pybamm-cspsur_ref)';
csnsur_rmse_table.pybamm = rms(csnsur_pybamm-csnsur_ref)';
cspave_rmse_table.pybamm = rms(cspave_pybamm-cspave_ref)';
csnave_rmse_table.pybamm = rms(csnave_pybamm-csnave_ref)';


%============= voltage ==============%
figure;
subplot(3,1,1)
plot(N_DIMS_LIST,V_rmse_table.pybamm,'--','Linewidth',2)
hold on
plot(N_DIMS_LIST,V_rmse_table.fvm,'-.','Linewidth',2)
hold on
plot(N_DIMS_LIST,V_rmse_table.fdm,':','Linewidth',2)

xlabel('Time [s]');
ylabel('Voltage [V]');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');
legend('Pybamm','Matlab-fvm','Matlab-fdm')

disp(V_rmse_table);


%============= average concentration ==============%

subplot(3,1,2)
plot(N_DIMS_LIST,cspave_rmse_table.pybamm,'--','Linewidth',2)
hold on
plot(N_DIMS_LIST,cspave_rmse_table.fvm,'-.', 'Linewidth',2)
hold on
plot(N_DIMS_LIST,cspave_rmse_table.fdm,':','Linewidth',2)

xlabel('Time [s]');
ylabel('Concentration [mol/m^3]');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');
legend('Pybamm','Matlab-fvm','Matlab-fdm')

disp(cspave_rmse_table);



subplot(3,1,3)
plot(N_DIMS_LIST,csnave_rmse_table.pybamm,'--','Linewidth',2)
hold on
plot(N_DIMS_LIST,csnave_rmse_table.fvm,'-.','Linewidth',2)
hold on
plot(N_DIMS_LIST,csnave_rmse_table.fdm,':','Linewidth',2)

xlabel('Time [s]');
ylabel('Concentration [mol/m^3]');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');
legend('Pybamm','Matlab-fvm','Matlab-fdm')

disp(csnave_rmse_table);

set(gcf,'unit','centimeters','position',[1,1.5,20,20])



%============= surface concentration ==============%
figure;
subplot(2,1,1)
plot(N_DIMS_LIST,cspsur_rmse_table.pybamm,'--','Linewidth',2)
hold on
plot(N_DIMS_LIST,cspsur_rmse_table.fvm,'-.','Linewidth',2)
hold on
plot(N_DIMS_LIST,cspsur_rmse_table.fdm,':','Linewidth',2)

xlabel('Time [s]');
ylabel('Concentration [mol/m^3]');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');
legend('Pybamm','Matlab-fvm','Matlab-fdm')

disp(cspsur_rmse_table);


subplot(2,1,2)
plot(N_DIMS_LIST,csnsur_rmse_table.pybamm,'--','Linewidth',2)
hold on
plot(N_DIMS_LIST,csnsur_rmse_table.fvm,'-.','Linewidth',2)
hold on
plot(N_DIMS_LIST,csnsur_rmse_table.fdm,':','Linewidth',2)

xlabel('Time [s]');
ylabel('Concentration [mol/m^3]');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');
legend('Pybamm','Matlab-fvm','Matlab-fdm')

disp(csnsur_rmse_table);

