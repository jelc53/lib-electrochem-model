function [De_n,De_s,De_p] = diffusivity_update(ce_n,ce_s,ce_p,T_core,param)

%Update Diffusivity in each region based on changing concentration
%Define molarities of each electrolyte region [mol/L] (aka mol/dm) (convert from mol/m3)

M_n = ce_n/1000; 
M_s = ce_s/1000;
M_p = ce_p/1000;

%Create vectors of core temperatures to match dimensions of Concentration vectors
T_c_n = repelem(T_core,param.Nx_n)';
T_c_s = repelem(T_core,param.Nx_s)';
T_c_p = repelem(T_core,param.Nx_p)';

% % From Tanim et al., 2014 (https://doi.org/10.1115/1.4028154)
% %NOTE: De_scaling factor included for sensitivity analysis
% De_n = 10.^(-(4.43 + 54./(T_c_n - (229+5*M_n)) + 0.22*M_n))*(1/1e4); %[m2/s]
% 
% De_s = 10.^(-(4.43 + 54./(T_c_s - (229+5*M_s)) + 0.22*M_s))*(1/1e4); %[m2/s]
% 
% De_p = 10.^(-(4.43 + 54./(T_c_p - (229+5*M_p)) + 0.22*M_p))*(1/1e4); %[m2/s]

% From Chen et al., 2020 (https://doi.org/10.1149/1945-7111/ab9050)
De_n = 8.794e-11*M_n.^2 - 3.972e-10*M_n + 4.862e-10; %[m2/s]

De_s = 8.794e-11*M_s.^2 - 3.972e-10*M_s + 4.862e-10; %[m2/s]

De_p = 8.794e-11*M_p.^2 - 3.972e-10*M_p + 4.862e-10; %[m2/s]

end