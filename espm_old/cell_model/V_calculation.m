function [V_cell, R_l] = V_calculation(ocp_p,ocp_n,eta_p,eta_n,phi_e,R_sei,I,SOC_cc, param)

% Include lumped resistance as function of SOC, using Coulomb-counted SOC
% Linear model Poly6 from R0_identification.m
%      f(x) = p1*x^6 + p2*x^5 + p3*x^4 + p4*x^3 + p5*x^2 + 
%                     p6*x + p7
% Coefficients (with 95% confidence bounds):
p1 =       1.251;  %(-0.168, 2.669)
p2 =      -4.205;  %(-8.467, 0.05748)
p3 =       5.629;  %(0.6152, 10.64)
p4 =      -3.824;  %(-6.741, -0.9083)
p5 =       1.392;  %(0.5237, 2.261)
p6 =     -0.2598;  %(-0.3822, -0.1374)
p7 =     0.04668;  %(0.04055, 0.05281)

% NOTE: Coefficients p_i, i=1...7 have been identified offline from HPPC at
% different SOC for the battery tested here (NMC, 4.85Ah) [Pozzato, Data in Brief, 2022]

R_l = p1.*SOC_cc.^6 + p2.*SOC_cc.^5 + p3.*SOC_cc.^4 + p4.*SOC_cc.^3 + ...
            p5.*SOC_cc.^2 + p6.*SOC_cc + p7;

% Cell voltage
% V_cell = ocp_p - ocp_n + eta_p - eta_n + phi_e - R_l.*I - R_sei*I;
% V_cell = ocp_p - ocp_n + eta_p - eta_n + phi_e;
V_cell = ocp_p - ocp_n + eta_p - eta_n;  %TODO: +phi_e

% Fix errant values to minimum voltage
% for i = 1:param.Nc
%     if V_cell(i)<2.50
%         V_cell(i) = 2.50;
%     else
%         V_cell(i) = V_cell(i);
%     end
%     if isreal(V_cell(i))
%         V_cell(i) = V_cell(i);
%     else
%         V_cell(i) = 2.50;
%     end      
% end

end