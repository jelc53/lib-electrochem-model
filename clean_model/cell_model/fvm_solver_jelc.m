function [t_out,u_out,te,xe,ie] = fvm_solver_jelc(x_initial, param, tspan)

t_out=tspan'; te=[]; xe=[]; ie=[];

%% Finite Volume implementation
u_out = x_initial(2:end-1);
u_i = u_out; u_outi = [];
t_steps = 1 / param.dt;

for t = 1:(param.t_duration / param.dt)
    % cathode
    [u_pi] = fvm_update_function(u_i(1:param.Nr-2), param, ...
        param.delta_xp, param.Dsp_ref, param.F, param.A, param.Lp, param.a_sp);
    
    % anode
    % u_out(param.Nr-1:end,t)
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