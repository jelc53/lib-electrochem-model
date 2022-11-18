function [f] = fvm_viscous_flux(r, dr, u_i, u_im1)
    f = (r^2 / dr) * (u_im1 - u_i);
end  % r^2? 