clear all
clc

all_matlab_fdm=load('matlab_results_fdm.mat');
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
plot(t_fdm,V_fdm,'--','Linewidth',2)

xlabel('Time [s]');
ylabel('Voltage [V]');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');
legend('Matlab-fdm','Pybamm')



rmseV=rms(V_fdm-V_pybamm)*1000;
disp(['RMSE of Voltage fdm=',num2str(rmseV),' mV'])


%=============average concentraion

subplot(3,1,2)
plot(t_fdm,cspave_fdm,'Linewidth',2)
hold on
plot(t_pybamm,cspave_pybamm,'--','Linewidth',2)
xlabel('Time [s]');
ylabel('Concentration [mol/m^3]');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');
legend('Matlab-fdm','Pybamm')



subplot(3,1,3)
plot(t_fdm,csnave_fdm,'Linewidth',2)
hold on
plot(t_pybamm,csnave_pybamm,'--','Linewidth',2)

xlabel('Time [s]');
ylabel('Concentration [mol/m^3]');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');
legend('Matlab-fdm','Pybamm')



rmse_avep=rms(cspave_fdm-cspave_pybamm);
rmse_aven=rms(csnave_fdm-csnave_pybamm);
disp(['RMSE of csp_ave fdm =',num2str(rmse_avep),' mol/m^3'])
disp(['RMSE of csn_ave fdm=',num2str(rmse_aven),' mol/m^3'])

set(gcf,'unit','centimeters','position',[1,1.5,20,20])







%=============Surface concentraion
figure;
subplot(2,1,1)
plot(t_fdm,cspsur_fdm,'Linewidth',2)
hold on
plot(t_pybamm,cspsur_pybamm,'--','Linewidth',2)
xlabel('Time [s]');
ylabel('Concentration [mol/m^3]');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');
legend('Matlab-fdm','Pybamm')


subplot(2,1,2)
plot(t_fdm,csnsur_fdm,'Linewidth',2)
hold on
plot(t_pybamm,csnsur_pybamm,'--','Linewidth',2)

xlabel('Time [s]');
ylabel('Concentration [mol/m^3]');
fontsize=14;
set(gca,'linewidth',1,'fontsize',fontsize,'fontname','Times New Roman');
legend('Matlab-fdm','Pybamm')


