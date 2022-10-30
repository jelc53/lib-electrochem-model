function [f] = viscous_flux(r, dr, u_i, u_im1)
    f = (r^2 / dr) * (u_i - u_im1);
end