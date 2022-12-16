function [Dsn, Dsp, kn, kp, Dsolv] = arrhenius_temp(param, T_core)

%Removed dependence on temperature

Dsn = param.Dsn_ref;
Dsp = param.Dsp_ref;
kn = param.kn_ref;
kp = param.kp_ref;
Dsolv = param.Dsolv_ref;

end