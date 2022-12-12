clear all
clc

all_matlab_fvm=load('matlab_results_fvm.mat');
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
    
    







t_pybamm=load('pybamm_t.txt');
V_pybamm=load('pybamm_V_cycle.txt');


cspave_pybamm=load('pybamm_cspave_cycle.txt');
csnave_pybamm=load('pybamm_csnave_cycle.txt');


cspsur_pybamm=load('pybamm_cspsur_cycle.txt');
csnsur_pybamm=load('pybamm_csnsur_cycle.txt');



figure;
subplot(3,1,1)
plot(t_pybamm,V_pybamm,'Linewidth',2)
hold on
plot(t_fvm,V_fvm,'--','Linewidth',2)

xlabel('Time [s]');
ylabel('Voltage [V]');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');
legend('Matlab-fvm','Pybamm')



rmseV=rms(V_fvm-V_pybamm)*1000;
disp(['RMSE of Voltage fvm=',num2str(rmseV),' mV'])


%=============average concentraion

subplot(3,1,2)
plot(t_fvm,cspave_fvm,'Linewidth',2)
hold on
plot(t_pybamm,cspave_pybamm,'--','Linewidth',2)
xlabel('Time [s]');
ylabel('Concentration [mol/m^3]');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');
legend('Matlab-fvm','Pybamm')



subplot(3,1,3)
plot(t_fvm,csnave_fvm,'Linewidth',2)
hold on
plot(t_pybamm,csnave_pybamm,'--','Linewidth',2)

xlabel('Time [s]');
ylabel('Concentration [mol/m^3]');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');
legend('Matlab-fvm','Pybamm')



rmse_avep=rms(cspave_fvm-cspave_pybamm);
rmse_aven=rms(csnave_fvm-csnave_pybamm);
disp(['RMSE of csp_ave fvm =',num2str(rmse_avep),' mol/m^3'])
disp(['RMSE of csn_ave fvm=',num2str(rmse_aven),' mol/m^3'])

set(gcf,'unit','centimeters','position',[1,1.5,20,20])







%=============Surface concentraion
figure;
subplot(2,1,1)
plot(t_fvm,cspsur_fvm,'Linewidth',2)
hold on
plot(t_pybamm,cspsur_pybamm,'--','Linewidth',2)
xlabel('Time [s]');
ylabel('Concentration [mol/m^3]');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');
legend('Matlab-fvm','Pybamm')


subplot(2,1,2)
plot(t_fvm,csnsur_fvm,'Linewidth',2)
hold on
plot(t_pybamm,csnsur_pybamm,'--','Linewidth',2)

xlabel('Time [s]');
ylabel('Concentration [mol/m^3]');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');
legend('Matlab-fvm','Pybamm')


