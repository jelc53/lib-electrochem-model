%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Evaluate stability & accuracy of model %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
warning off
clc

% Voltage
ref = 'nr100_cr1'
rms_voltage = dictionary()
results_voltage = table() 
for nr = [100, 50, 10]
    for cr = [1]
        nr_str = string(nr);
        cr_str = string(cr);
        identifier = 'nr' + nr_str + '_cr' + cr_str
        filename = 'output/fdm_ode15s_' + identifier + '_voltage.dat';
        out_vec = transpose(csvread(filename));
        results_voltage = addvars(results_voltage, out_vec, 'NewVariableNames', identifier);
        
        % rms calc
        metric = (1/length(out_vec)) * sqrt(sum((results_voltage.(ref) - out_vec).^2))
        rms_voltage(identifier) = metric
    end
end

% Surface concentration
ref = 'nr100_cr1'
rms_csurf = dictionary()
results_csurf = table() 
for nr = [100, 50, 10]
    for cr = [1]
        nr_str = string(nr);
        cr_str = string(cr);
        identifier = 'nr' + nr_str + '_cr' + cr_str
        filename = 'output/fdm_ode15s_' + identifier + '_csurf.dat';
        out_vec = transpose(csvread(filename));
        results_csurf = addvars(results_csurf, out_vec, 'NewVariableNames', identifier);

        % rms calc
        metric = (1/length(out_vec)) * sqrt(sum((results_csurf.(ref) - out_vec).^2))
        rms_csurf(identifier) = metric
    end
end
