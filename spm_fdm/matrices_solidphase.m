function [A_sd, B_sd_n, B_sd_p] = matrices_solidphase(param)
%%%%%% ODE matrix formulation for Solid Phase of a cell %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% ODE: State-space formulation
A_mat = zeros(param.Nr-1);
B_mat = zeros(param.Nr-1,1);
for k = 1:(param.Nr-1)
    if k == 1
        A_mat(k,k) = -2;
        A_mat(k,k+1) = 2;
        B_mat(k,1) = 0;
    elseif k == (param.Nr-1)
        A_mat(k,k-1) = 2;
        A_mat(k,k) = -2;
        B_mat(k,1) = 2*(1+1/k);
    else
        A_mat(k,k-1) = 1*(1-1/k);
        A_mat(k,k) = -2;
        A_mat(k,k+1) = 1*(1+1/k);
        B_mat(k,1) = 0;
    end
end

A_sd = sparse(A_mat); 

% Creating a column vector with B vector for all cells
beta_n = 1/(param.F*param.A*param.Ln*param.a_sn*param.delta_xn);
beta_p = 1/(param.F*param.A*param.Lp*param.a_sp*param.delta_xp);

B_sd_n = sparse(-B_mat*beta_n);
B_sd_p = sparse(B_mat*beta_p);

end