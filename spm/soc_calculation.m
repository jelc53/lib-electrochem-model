function [soc_bulk_n,soc_bulk_p] = soc_calculation(cs_n,cs_p,param)

soc_int_n = [param.delta_xn:param.delta_xn:param.delta_xn*(param.Nr-1)]';
soc_int_p = [param.delta_xp:param.delta_xp:param.delta_xp*(param.Nr-1)]';

% for j = 1:length(param.t_data)
for j = 1:length(cs_n) %NOTE: Changed this from param.t_data for ODE event function handling
for i = 1:param.Nc
shift_R = (i-1)*(param.Nr-1);
xbulk_n(i,j) = trapz(soc_int_n,soc_int_n.^2.*cs_n(shift_R+1:shift_R+...
    (param.Nr-1),j))*3/(param.c_n_max*param.Rs_n^3);
xbulk_p(i,j) = trapz(soc_int_p,soc_int_p.^2.*cs_p(shift_R+1:shift_R+...
    (param.Nr-1),j))*3/(param.c_p_max*param.Rs_p^3);

soc_bulk_n(i,j) = (xbulk_n(i,j) - param.theta0_n)/(param.theta100_n - param.theta0_n);
soc_bulk_p(i,j) = (param.theta0_p - xbulk_p(i,j))/(param.theta0_p - param.theta100_p);

if soc_bulk_n(i,j)>1
    soc_bulk_n(i,j) = 1;
end
if soc_bulk_p(i,j)<0
    soc_bulk_p(i,j) = 0;
end
    if soc_bulk_n(i,j)>1
        soc_bulk_n(i,j) = 1;
    elseif soc_bulk_n(i,j)<0
        soc_bulk_n(i,j) = 0;
    else
        soc_bulk_n(i,j) = soc_bulk_n(i,j);           
    end
    if soc_bulk_p(i,j)>1
        soc_bulk_p(i,j) = 1;
    elseif soc_bulk_p(i,j)<0
        soc_bulk_p(i,j) = 0;
    else
        soc_bulk_p(i,j) = soc_bulk_p(i,j);           
    end

end
end

end
