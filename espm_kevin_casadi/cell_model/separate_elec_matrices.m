%% Separate compiled coefficient matrices - Finite Volume Coefficient Scheme
param.A1_en_d = param.A_en_d(1:param.Nx_n*param.Nc,1:param.Nx_n*param.Nc);
param.A2_en_d = param.A_en_d(param.Nx_n*param.Nc+1:2*param.Nx_n*param.Nc,...
    1:param.Nx_n*param.Nc);

param.A1_es_d = param.A_es_d(1:param.Nx_s*param.Nc,1:param.Nx_s*param.Nc);
param.A2_es_d = param.A_es_d(param.Nx_s*param.Nc+1:2*param.Nx_s*param.Nc,...
    1:param.Nx_s*param.Nc);

param.A1_ep_d = param.A_ep_d(1:param.Nx_p*param.Nc,1:param.Nx_p*param.Nc);
param.A2_ep_d = param.A_ep_d(param.Nx_p*param.Nc+1:2*param.Nx_p*param.Nc,...
    1:param.Nx_p*param.Nc);

param.B_en_d = param.B_e_d(1:param.Nx_n*param.Nc,1);
param.B_es_d = param.B_e_d(param.Nx_n*param.Nc+1:(param.Nx_n+param.Nx_s)*param.Nc,1);
param.B_ep_d = param.B_e_d((param.Nx_n+param.Nx_s)*param.Nc+1:end,1);
