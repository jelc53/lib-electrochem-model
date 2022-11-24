function [f] = fvm_viscous_flux(r, dr, u_i, u_im1)
    % convert to concentration
    c_i = u_i / (r+dr/2)^2;
    c_im1 = u_im1 / (r-dr/2)^2;
    
    % compute gradient
    f = (1 / dr) * (c_i - c_im1);
end 