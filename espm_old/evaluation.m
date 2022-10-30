%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Evaluate stability & accuracy of model %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
warning off
clc

model = 'spm_fdm';
ref_model = 'espm_fdm';

%% Reference solutions
% load data from previous runs
ref_iapp24 = load(fullfile('output',append(ref_model,'_nr1000_cr1_time3600_cyc0.mat')));
ref_iapp48 = load(fullfile('output',append(ref_model,'_nr1000_cr2_time1000_cyc0.mat')));
ref_iapp96 = load(fullfile('output',append(ref_model,'_nr1000_cr4_time300_cyc0.mat')));
ref_data = [ref_iapp24 ref_iapp48 ref_iapp96];

ref_table = table;
ref_cs_avg_n = [];
ref_cs_avg_p = [];
for i = 1:length(ref_data)
    ref = ref_data(i);
    Nr = ref.all_data.node_Nr;
    volume_n = (4/3)*pi*(ref.all_data.par.Rs_n^3);
    volume_p = (4/3)*pi*(ref.all_data.par.Rs_p^3);
    delta_xn = ref.all_data.par.delta_xn;
    delta_xp = ref.all_data.par.delta_xp;
    vol_xn = (4/3)*pi*((delta_xn * (1:(Nr-1))).^3);
    vol_xp = (4/3)*pi*((delta_xp * (1:(Nr-1))).^3);
    weights_n = [vol_xn(1), diff(vol_xn)] ./ volume_p;
    weights_p = [vol_xn(1), diff(vol_xp)] ./ volume_p;
    cs_avg_n = dot(weights_n, ref.all_data.csn(:,end));
    cs_avg_p = dot(weights_p, ref.all_data.csp(:,end));
    
    ref_cs_avg_n = [ref_cs_avg_n; cs_avg_n];
    ref_cs_avg_p = [ref_cs_avg_p; cs_avg_p];

    ref_table.(string(i)) = [
        ref.all_data.V(end); 
        ref.all_data.csn(end,end);
        ref.all_data.csp(end,end);
        dot(weights_n, ref.all_data.csn(:,end));
        dot(weights_p, ref.all_data.csp(:,end));
    ];
end

ref_table.Properties.RowNames = {'voltage';'cathode_surf_con';'anode_surf_con';'cathode_avg_con';'anode_avg_con'};
ref_table.Properties.VariableNames = {'iapp24';'iapp48';'iapp96'};
writetable(ref_table, fullfile('ref_table.dat'))


%% Experiments to evaluate accuracy
results_table = table;
c_rates = [1 2 4];
nr_dims = [5 10 100 1000];
cr_duration = [3600 1000 300];
variable_names = {'measure','nr','c_rate','t_duration','last_ref','last_val','rmse'};

for i = 1:length(c_rates)
    for j = 1:length(nr_dims)
        % load experiment data
        filename = string(model ) + '_nr' + string(nr_dims(j)) + '_cr' + string(c_rates(i)) + '_time' + string(cr_duration(i)) + '_cyc0.mat';
        exp_data = load(fullfile('output',filename));

        % voltage
        ref_vals = ref_data(i).all_data.V;
        exp_vals = exp_data.all_data.V;
        new_row = {
            'voltage',...
            nr_dims(j),...
            c_rates(i),...
            cr_duration(i),...
            ref_vals(end),...
            exp_vals(end),...
            rmse(ref_vals,exp_vals)
        };
        results_table = [results_table; new_row];

        % surface concentration

        % average concentration
    end
end
results_table.Properties.VariableNames = variable_names;
writetable(results_table, fullfile('results_table.dat'))

%% Plot V_RMS, Csurf_RMS, Cavg_RMS
% Voltage
figure; hold on;  grid on;
a = [];
for i = 1:length(c_rates)
    idx = strcmp(results_table.measure{i},'voltage') & results_table.c_rate == c_rates(i);
    a = [a, semilogx(results_table.nr(idx), results_table.rmse(idx), 'linewidth', 1)];
end
M1 = "c-rate = 1"; M2 = "c-rate = 2"; M3 = "c-rate = 4";
legend(a, [M1, M2, M3]);
xlabel('Spatial resolution (Nr)','Fontsize',12,'interpreter','latex')
ylabel('V_RMSE','Fontsize',12,'interpreter','latex')
set(gca,'Fontsize',12, 'XScale', 'log');
saveas(gcf,'voltage_rmse_plot.png')