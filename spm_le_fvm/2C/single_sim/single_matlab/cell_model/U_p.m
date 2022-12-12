function y = U_p(x)

global err_Up err_theta_surf_p

% [chen] Development of Experimental Techniques for Parameterization of Multiscale Lithium-ion Battery Models

% Correction error
if length(err_Up) > 1
    err = interp1(err_theta_surf_p,err_Up,x);   
    if isnan(err)
        err = 0;
    end
else
    err = 0;
end
     
y = -0.8090*x + 4.4875 - 0.0428*tanh(18.5138*(x-0.5542))...
    -17.7326*tanh(15.7890*(x-0.3117))+17.5842*tanh(15.9308*(x-0.3120))...
    +err;
end