function [u_i] = fvm_update_function(u_out, param, dr, D_sj, F, A, L_j, a_sj)
    u_i = zeros(1,param.Nr-2);

    for k = 1:(param.Nr-2)
        r = dr*k - (dr/2);

        % compute flux plus and minus
        if k == 1  % centre of particle
            f_plus = fvm_viscous_flux(r, dr, u_out(k+1), u_out(k));
            f_minus = fvm_viscous_flux(r, dr, u_out(k), u_out(k));

        elseif k == param.Nr-2  % surface of particle
            f_plus = r^2 * param.I_data(1) / (D_sj*F*A*L_j*a_sj);  % constant current
            f_minus = fvm_viscous_flux(r, dr, u_out(k), u_out(k-1));

        else
            f_plus = fvm_viscous_flux(r, dr, u_out(k+1), u_out(k));
            f_minus = fvm_viscous_flux(r, dr, u_out(k), u_out(k-1));

        end
        
        % conservative update schematic (forward euler)
        u_i(k) = u_out(k) - D_sj * (param.dt / dr) * (f_plus - f_minus);

    end

end