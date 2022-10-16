%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Evaluate stability & accuracy of model %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
warning off
clc

% Voltage
results_voltage = table() 
for nr = [10, 50, 100]
    for cr = [1]
        nr_str = string(nr);
        cr_str = string(cr);
        identifier = 'nr' + nr_str + '_cr' + cr_str
        filename = 'output/fdm_ode15s_' + identifier + '_voltage.dat';
        out_vec = transpose(csvread(filename));
        results_voltage = addvars(results_voltage, out_vec, 'NewVariableNames', identifier);
    end
end

% Surface concentration
results_csurf = table() 
for nr = [10, 50, 100]
    for cr = [1]
        nr_str = string(nr);
        cr_str = string(cr);
        identifier = 'nr' + nr_str + '_cr' + cr_str
        filename = 'output/fdm_ode15s_' + identifier + '_csurf.dat';
        out_vec = transpose(csvread(filename));
        results_csurf = addvars(results_csurf, out_vec, 'NewVariableNames', identifier);
    end
end
