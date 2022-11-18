function [t_out,cell_out,te,xe,ie] = fvm_solver_jelc(x_initial, param, tspan)

    t_out=tspan'; te=[]; xe=[]; ie=[];
    
    %% Finite Volume implementation
    cell_out = x_initial(2:end-1);
    rp_sq_vec = generate_radius_vec(param.delta_xp, param.Nr).^2;
    rn_sq_vec = generate_radius_vec(param.delta_xn, param.Nr).^2;
    u_out = cell_out.*[rp_sq_vec rn_sq_vec]';
    
    u_i = u_out; u_outi = [];
    t_steps = 1 / param.dt;

    for t = 1:(param.t_duration / param.dt)
        % cathode
        [u_pi] = fvm_update_function(u_i(1:param.Nr-2), param, ...
            param.delta_xp, param.Dsp_ref, param.F, param.A, param.Lp, param.a_sp);
        
        % anode
        [u_ni] = fvm_update_function(u_i(param.Nr-1:end), param, ...
            param.delta_xn, param.Dsn_ref, param.F, param.A, param.Ln, param.a_sn);
        
        % combine and update 
        u_i = [u_pi, u_ni];
        u_outi = [u_outi u_i'];
        
        if mod(t,t_steps) == 0
            u_out = [u_out mean(u_outi,2)]; 
            u_outi = [];
        end

    end
    cell_out = u_out .* [1./rp_sq_vec 1./rp_sq_vec]';
end
