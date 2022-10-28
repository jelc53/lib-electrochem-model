function [eps_el_n, R_sei] = sei_lam_effect(L_film,param,aina_n,af_n)

% Modified electrolyte porosity
% eps_el_n = 1-param.epsilon_n*(1+3*L_film./param.Rs_n)-param.eps_filler_n + (aina_n-af_n)*param.Rs_n/3;
eps_el_n = 1-param.epsilon_n-param.eps_filler_n;
% Deactivate SEI layer resistance
R_sei = 1/((param.a_sn+af_n-aina_n)*param.A*param.Ln*param.kappa_sei) * 0;
% 
% % SEI layer resistance
% R_sei = 1/((param.a_sn+af_n-aina_n)*param.A*param.Ln*param.kappa_sei);

end