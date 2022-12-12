function [Alam_mat_n,Alam_mat_p,Blam_mat_n,Blam_mat_p] = matrices_lam(param)

% ODE matrix formulation of LAM for cathode and anode %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% ODE: State-space formulation
Alam_mat_n = -ones(param.Nc)*param.beta_lam_n;
Alam_mat_p = -ones(param.Nc)*param.beta_lam_p;

b_lam_n = param.beta_lam_n*[1 param.k_lam_n];
b_lam_p = param.beta_lam_p*[1 param.k_lam_p];
Blam_mat_n = repmat(b_lam_n*param.a_sn,param.Nc,1);
Blam_mat_p = repmat(b_lam_p*param.a_sp,param.Nc,1);

end