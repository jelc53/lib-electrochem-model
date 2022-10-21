function [r_i, u_i] = finite_volume_update(r_out, u_out, param, dr, D_sj, F, A, L_j, a_sj)
    u_i = zeros(1,param.Nr-2);

    for k = 1:(param.Nr-2)
        r = dr*k - (dr/2);

        % compute flux plus and minus
        if k == 1
            f_plus = viscous_flux(r, dr, u_out(k+1), u_out(k));
            f_minus = viscous_flux(r, dr/2, u_out(k), r_out(k));

        elseif k == param.Nr-2
            f_plus = viscous_flux(r, dr/2, r_out(k+1), u_out(k));
            f_minus = viscous_flux(r, dr, u_out(k), u_out(k-1));

        else
            f_plus = viscous_flux(r, dr, u_out(k+1), u_out(k));
            f_minus = viscous_flux(r, dr, u_out(k), u_out(k-1));

        end
        
        % conservative update schematic
        u_i(k) = u_out(k) - (param.dt / dr) * (f_plus - f_minus);

    end
%     u_out = [u_out; u_i];

    r_i = movmean(u_i, 2);
    r_i(1) = r_out(1);  % no change for radius = 0
    r_i_surf = param.I_data(1) / (D_sj*F*A*L_j*a_sj);  % constant current
    r_i = [r_i, r_i_surf];
%     r_out = [r_out; r_i];
end