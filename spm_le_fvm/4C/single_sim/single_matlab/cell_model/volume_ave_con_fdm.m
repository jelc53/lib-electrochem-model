function [csp_ave,csn_ave] = volume_ave_con_fdm(param,cs_p,cs_n)

r_n = [param.delta_xn:param.delta_xn:param.delta_xn*(param.Nr-1)]';
r_p = [param.delta_xp:param.delta_xp:param.delta_xp*(param.Nr-1)]';

csn_ave = trapz(r_n,r_n.^2.*cs_n(1:(param.Nr-1),:))*3/(param.Rs_n^3);
csp_ave = trapz(r_p,r_p.^2.*cs_p(1:(param.Nr-1),:))*3/(param.Rs_p^3);

end

