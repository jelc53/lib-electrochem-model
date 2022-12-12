function [dus_dt] = fvm_ode(t_in, u_in, param)

    % anode density
    [dusn_dt] = fvm_update_function(1, ceil(t_in), u_in(1:param.Nr-1), param, ...
        -1, param.delta_xn, param.Dsn_ref, param.F, param.A, param.Ln, param.a_sn);

    % cathode density
    [dusp_dt] = fvm_update_function(1, ceil(t_in), u_in(param.Nr:end), param, ...
        1, param.delta_xp, param.Dsp_ref, param.F, param.A, param.Lp, param.a_sp);
    
    % combine and update 
    dus_dt = [dusn_dt dusp_dt]';

end