function [cs_initial, csn0, csp0] = conc_initial_sd(SOC_cells, param)

csn0 = (SOC_cells*(param.theta100_n - param.theta0_n) + param.theta0_n).*param.c_n_max ;
csp0 = (param.theta0_p - SOC_cells*(param.theta0_p - param.theta100_p)).*param.c_p_max; 

csn_initial = zeros(param.Nc*(param.Nr-1),1);
csp_initial = zeros(param.Nc*(param.Nr-1),1);

% Initialize concentration states in both electrodes
for i=1:param.Nc
    csn_initial((i-1)*(param.Nr-1)+1:i*(param.Nr-1),1) = repmat(csn0(i),(param.Nr-1),1);
    csp_initial((i-1)*(param.Nr-1)+1:i*(param.Nr-1),1) = repmat(csp0(i),(param.Nr-1),1);
end

cs_initial = [csn_initial;csp_initial];

end