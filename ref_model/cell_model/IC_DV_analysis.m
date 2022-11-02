function [DV,IC] = IC_DV_analysis(cutoff_freq,Ts,V_cell,disch_capacity)
% ICA and DVA
%
% coff_freq [Hz]
% V_cell [V]
% disch_capacity [Ah]

% Filter 
s = tf('s');
f = cutoff_freq;                                                            % [Hz]
damp = 1;                                                                   % [-] damping

% Bode filter
fdt_tmp = (2*pi*f)^2/(s^2 + 2*pi*f*s*damp + (2*pi*f)^2);

% DV
DV_raw = diff(V_cell)./diff(disch_capacity);
DV_raw(find(isinf(DV_raw))) = 10;

DV = my_filter(DV_raw,'lp',f,1/(2*damp),Ts,0,0);

% IC
IC = 1./DV;

end