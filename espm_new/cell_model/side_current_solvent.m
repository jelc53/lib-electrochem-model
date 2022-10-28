function [i_s] = side_current_solvent(R_sei, T_core, phi_sn, param, input_crt, theta_surf_n,Csolv)
import casadi.*
i_s=SX.sym('i_s',param.Nc);
kf = 2*param.k_f*exp(-param.Ea_kf*(1./T_core - 1/param.Tref)./param.Rg)*param.c_n_max^2.*...
    (theta_surf_n.^2).*exp(param.beta_val*param.F*param.U_s./(param.Rg*T_core));

for i = 1:param.Nc
    index = (i-1)*param.Nsei+1;

    i_s(i,1) = -param.F*kf(i)*Csolv(index).*exp(-param.beta_val*param.F./(param.Rg*T_core(i)).*...
        (phi_sn(i)-R_sei(i)*input_crt));
end

end