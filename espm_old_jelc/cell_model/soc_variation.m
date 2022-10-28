function [SOC_cells] = soc_variation(param)   

% Randomize the variation in SOC between cells in a series module

SOC_cells = param.SOC_ref*ones(param.Nc,1) + (param.SOC_deviation)*randn(param.Nc,1);
SOC_cells = soc_bound_module(SOC_cells);

end