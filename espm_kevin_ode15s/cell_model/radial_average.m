function [wj_vec, csj_avg] = radial_average(csj,Nr,Rs_j,delta_xj)

    vol_total = (4/3)*pi*(Rs_j^3);
    vol_radii_vec = (4/3)*pi*((delta_xj * (1:(Nr-1))).^3);
    wj_vec = [vol_radii_vec(1), diff(vol_radii_vec)] ./ vol_total;
    %csj_avg = dot(wj_vec, csj(:,end));
    csj_avg = wj_vec * csj;
  