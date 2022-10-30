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
ref_csn_avg = [];
ref_csp_avg = [];
for i = 1:length(ref_data)
    ref = ref_data(i);
    Nr = ref.all_data.node_Nr;
    Rs_n = ref.all_data.par.Rs_n; Rs_p = ref.all_data.par.Rs_p; 
    delta_xn = ref.all_data.par.delta_xn; delta_xp = ref.all_data.par.delta_xp; 
    [wn_vec, csn_avg] = radial_average(ref.all_data.csn(:,end),Nr,Rs_n,delta_xn);
    [wp_vec, csp_avg] = radial_average(ref.all_data.csp(:,end),Nr,Rs_p,delta_xp);

    ref_csn_avg = [ref_csn_avg; csn_avg];
    ref_csp_avg = [ref_csp_avg; csp_avg];

    ref_table.(string(i)) = [
        ref.all_data.V(end); 
        ref.all_data.csn(end,end);
        ref.all_data.csp(end,end);
        dot(wn_vec, ref.all_data.csn(:,end));
        dot(wp_vec, ref.all_data.csp(:,end));
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
        
        Nr = exp_data.all_data.node_Nr;
        ref_Nr = ref_data(i).all_data.node_Nr;
        Rs_n = exp_data.all_data.par.Rs_n; Rs_p = exp_data.all_data.par.Rs_p;
        delta_xn = exp_data.all_data.par.delta_xn;  delta_xp = exp_data.all_data.par.delta_xp;

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

        % surface concentration (anode)
        ref_vals = ref_data(i).all_data.csn(end,:);
        exp_vals = exp_data.all_data.csn(end,:);
        new_row = {
            'csn_surf',...
            nr_dims(j),...
            c_rates(i),...
            cr_duration(i),...
            ref_vals(end),...
            exp_vals(end),...
            rmse(ref_vals,exp_vals)
        };
        results_table = [results_table; new_row];

        % surface concentration (cathode)
        ref_vals = ref_data(i).all_data.csp(end,:);
        exp_vals = exp_data.all_data.csp(end,:);
        new_row = {
            'csp_surf',...
            nr_dims(j),...
            c_rates(i),...
            cr_duration(i),...
            ref_vals(end),...
            exp_vals(end),...
            rmse(ref_vals,exp_vals)
        };
        results_table = [results_table; new_row];

        % average concentration (anode)
        [ref_wn_vec, ref_vals] = radial_average(ref_data(i).all_data.csn,ref_Nr,Rs_n,delta_xn);
        [wn_vec, exp_vals] = radial_average(exp_data.all_data.csn,Nr,Rs_n,delta_xn);
        new_row = {
            'csn_avg',...
            nr_dims(j),...
            c_rates(i),...
            cr_duration(i),...
            ref_vals(end),...
            exp_vals(end),...
            rmse(ref_vals,exp_vals)
        };
        results_table = [results_table; new_row];

        % average concentration (anode)
        [ref_wp_vec, ref_vals] = radial_average(ref_data(i).all_data.csp,ref_Nr,Rs_p,delta_xp);
        [wp_vec, exp_vals] = radial_average(exp_data.all_data.csp,Nr,Rs_p,delta_xp);
        new_row = {
            'csp_avg',...
            nr_dims(j),...
            c_rates(i),...
            cr_duration(i),...
            ref_vals(end),...
            exp_vals(end),...
            rmse(ref_vals,exp_vals)
        };
        results_table = [results_table; new_row];
    end
end
results_table.Properties.VariableNames = variable_names;
writetable(results_table, fullfile('results_table.dat'))

%% Plot V_RMS, Csurf_RMS, Cavg_RMS
% voltage
figure; hold on;  grid on;
a = [];
for i = 1:length(c_rates)
    idx = strcmp(results_table.measure,'voltage') & results_table.c_rate == c_rates(i);
    a = [a, semilogx(results_table.nr(idx), results_table.rmse(idx), 'linewidth', 1)];
end
M1 = "c-rate = 1"; M2 = "c-rate = 2"; M3 = "c-rate = 4";
legend(a, [M1, M2, M3]);
xlabel('Spatial resolution (Nr)','Fontsize',12,'interpreter','latex')
ylabel('Voltage RMSE','Fontsize',12,'interpreter','latex')
set(gca,'Fontsize',12, 'XScale', 'log');
saveas(gcf,'voltage_rmse_plot.png')

% surface concentration (anode)
figure; hold on;  grid on;
a = [];
for i = 1:length(c_rates)
    idx = strcmp(results_table.measure,'csn_surf') & results_table.c_rate == c_rates(i);
    a = [a, semilogx(results_table.nr(idx), results_table.rmse(idx), 'linewidth', 1)];
end
M1 = "c-rate = 1"; M2 = "c-rate = 2"; M3 = "c-rate = 4";
legend(a, [M1, M2, M3]);
xlabel('Spatial resolution (Nr)','Fontsize',12,'interpreter','latex')
ylabel('Anode Csurf RMSE','Fontsize',12,'interpreter','latex')
set(gca,'Fontsize',12, 'XScale', 'log');
saveas(gcf,'cn_surf_rmse_plot.png')

% surface concentration (cathode)
figure; hold on;  grid on;
a = [];
for i = 1:length(c_rates)
    idx = strcmp(results_table.measure,'csp_surf') & results_table.c_rate == c_rates(i);
    a = [a, semilogx(results_table.nr(idx), results_table.rmse(idx), 'linewidth', 1)];
end
M1 = "c-rate = 1"; M2 = "c-rate = 2"; M3 = "c-rate = 4";
legend(a, [M1, M2, M3]);
xlabel('Spatial resolution (Nr)','Fontsize',12,'interpreter','latex')
ylabel('Cathode Csurf RMSE','Fontsize',12,'interpreter','latex')
set(gca,'Fontsize',12, 'XScale', 'log');
saveas(gcf,'cp_surf_rmse_plot.png')

% average concentration (anode)
figure; hold on;  grid on;
a = [];
for i = 1:length(c_rates)
    idx = strcmp(results_table.measure,'csn_avg') & results_table.c_rate == c_rates(i);
    a = [a, semilogx(results_table.nr(idx), results_table.rmse(idx), 'linewidth', 1)];
end
M1 = "c-rate = 1"; M2 = "c-rate = 2"; M3 = "c-rate = 4";
legend(a, [M1, M2, M3]);
xlabel('Spatial resolution (Nr)','Fontsize',12,'interpreter','latex')
ylabel('Anode Cavg RMSE','Fontsize',12,'interpreter','latex')
set(gca,'Fontsize',12, 'XScale', 'log');
saveas(gcf,'cn_avg_rmse_plot.png')

% average concentration (cathode)
figure; hold on;  grid on;
a = [];
for i = 1:length(c_rates)
    idx = strcmp(results_table.measure,'csp_avg') & results_table.c_rate == c_rates(i);
    a = [a, semilogx(results_table.nr(idx), results_table.rmse(idx), 'linewidth', 1)];
end
M1 = "c-rate = 1"; M2 = "c-rate = 2"; M3 = "c-rate = 4";
legend(a, [M1, M2, M3]);
xlabel('Spatial resolution (Nr)','Fontsize',12,'interpreter','latex')
ylabel('Cathode Cavg RMSE','Fontsize',12,'interpreter','latex')
set(gca,'Fontsize',12, 'XScale', 'log');
saveas(gcf,'cp_avg_rmse_plot.png')