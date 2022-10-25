function [r_i, u_i] = finite_volume_update(r_out, u_out, param, dr, D_sj, F, A, L_j, a_sj)
    u_i = zeros(1,param.Nr-2);

    for k = 1:(param.Nr-2)
        r = dr*k - (dr/2);

        % compute flux plus and minus
        if k == 1  % centre of particle
            f_plus = viscous_flux(r, dr, u_out(k+1), u_out(k));
            f_minus = viscous_flux(r, dr, u_out(k), u_out(k));

        elseif k == param.Nr-2  % surface of particle
            f_plus = param.I_data(1) / (D_sj*F*A*L_j*a_sj);  % constant current
            f_minus = viscous_flux(r, dr, u_out(k), u_out(k-1));

        else
            f_plus = viscous_flux(r, dr, u_out(k+1), u_out(k));
            f_minus = viscous_flux(r, dr, u_out(k), u_out(k-1));

        end
        
        % conservative update schematic
        u_i(k) = u_out(k) - (D_sj*4*pi) * (param.dt / dr) * (f_plus - f_minus);

    end
        r_i = [movmean(u_i, 2), u_i(end)];
%     u_out = [u_out; u_i];
%     r_out = [r_out; r_i];
end