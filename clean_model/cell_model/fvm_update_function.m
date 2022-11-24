function [u_i] = fvm_update_function(mode, t, u_out, param, sgn, dr, D_sj, F, A, L_j, a_sj)
    u_i = zeros(1,param.Nr-1);

    for k = 1:(param.Nr-1)

        % cell radii 
        r = dr*k - (dr/2);
        r_plus = r+dr/2;
        r_minus = r-dr/2;

        % compute spherical flux
        if k == 1  % centre of particle
            f_plus = (r_plus)^2 * fvm_viscous_flux(r_plus, dr, u_out(k+1), u_out(k));
            f_minus = 0; %(r_minus)^2 * fvm_viscous_flux(r_minus, dr, u_out(k), u_out(k));

        elseif k == param.Nr-1  % surface of particle
            f_plus = (r_plus)^2 * sgn*param.I_data(t+1) / (D_sj*F*A*L_j*a_sj);  % constant current
            f_minus = (r_minus)^2 * fvm_viscous_flux(r_minus, dr, u_out(k), u_out(k-1));

        else
            f_plus = (r_plus)^2 * fvm_viscous_flux(r_plus, dr, u_out(k+1), u_out(k));
            f_minus = (r_minus)^2 * fvm_viscous_flux(r_minus, dr, u_out(k), u_out(k-1));

        end
        
        % conservative update
        if mode == 0
            % density update (forward euler)
            u_i(k) = u_out(k) + D_sj * (param.dt / dr) * (f_plus - f_minus);
        else
            % ode du_dt update (ode15s)
            u_i(k) = D_sj * (1 / dr) * (f_plus - f_minus);
        end

    end

end