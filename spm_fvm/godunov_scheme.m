function [uj_np1] godunov_scheme(r, dr, u_i, u_im1)
    if u_i > u_im1
        % do something
        min(1,0);
    else
        % do something else
        max(1,0);
    end
    uj_np1 = 0;
end