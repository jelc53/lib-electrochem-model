function [f] = fvm_viscous_flux(r, dr, u_i, u_im1)
    c_i = u_i / (r+dr/2)^2;
    c_im1 = u_im1 / (r-dr/2)^2;
    
    f = (1 / dr) * (c_i - c_im1);
end 