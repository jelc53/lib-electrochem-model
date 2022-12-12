function [r_vec] = generate_radius_vec(dr, Nr)
    for k = 1:(Nr-1)
        r_vec(k) = dr*k - (dr/2);
    end
end