%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Evaluate stability (mass conservation) of model %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
warning off
clc

model = 'spm_fdm';

%% Populate results tables
c_rates = [1 2 4];
nr_dims = [5 10 100];
t_duration = [3600 1000 300];
results = cell(length(c_rates));

for i = 1:length(c_rates)
    n_steps = t_duration(i)*2*10;
    results_table = table();  %'Size',[n_steps,6]
    
    for j = 1:length(nr_dims)
        % load experiment data
        filename = string(model) + '_nr' + string(nr_dims(j)) + '_cr' + string(c_rates(i)) + '_time' + string(t_duration(i)) + '_cyc10.mat';
        exp_data = load(fullfile('output',filename));

        % compute radial averages
        Nr = exp_data.all_data.node_Nr;
        Rs_n = exp_data.all_data.par.Rs_n; Rs_p = exp_data.all_data.par.Rs_p;
        delta_xn = exp_data.all_data.par.delta_xn;  delta_xp = exp_data.all_data.par.delta_xp;
        [wn_vec, csn_avg] = radial_average(exp_data.all_data.csn,Nr,Rs_n,delta_xn);
        [wp_vec, csp_avg] = radial_average(exp_data.all_data.csp,Nr,Rs_p,delta_xp);
        
        % append to results table
        col_name1 = 'csn_nr'+string(nr_dims(j));
        col_name2 = 'csp_nr'+string(nr_dims(j)); 
        results_table.(col_name1) = csn_avg';
        results_table.(col_name2) = csp_avg';
    
    end
%     results_table.Properties.VariableNames = {'csn_nr5', 'csp_nr5', 'csn_nr10', 'csp_nr10', 'csn_nr100', 'csp_nr100'};
    writetable(results_table, fullfile('cr'+string(c_rates(i))+'_cavg_cyc10_table.dat'));
    results{i} = results_table;
end

%% Plot Cavg over 10 cycles for each c-rate
electrode = ['n', 'p'];
for e = electrode
    for i = 1:length(c_rates)
        % negative electrode 
        figure('position', [0 0 600 300]); hold on;  grid on; a = [];
        results_table = results{i};
        for j = 1:length(nr_dims)
            col_name = 'cs' + string(e) + '_nr' + string(nr_dims(j));
            x = 1:(t_duration(i)*2*10);    
            y = results_table.(col_name);
            a = [a, plot(x, y, 'linewidth', 1)];
        end
        col_name = 'cs' + string(e) + '_nr5';
        yline(max(results_table.(col_name)(1:(2*t_duration(i)))), 'r');
        yline(min(results_table.(col_name)(1:(2*t_duration(i)))), 'r');
        M1 = "Nr = 5"; M2 = "Nr = 10"; M3 = "Nr = 100";
        legend(a, [M1, M2, M3]);
        xlabel('time step (sec)','Fontsize',12,'interpreter','latex')
        ylabel('cs'+string(e)+' average @ c-rate of '+string(c_rates(i)),'Fontsize',12,'interpreter','latex')
        set(gca,'Fontsize',12);
        filename = strcat('stability', '_cs', string(e), '_cr', string(c_rates(i)), '_plot.png');
        saveas(gcf,fullfile('output', filename)) 
    end
end


%% TODO Plot V_cell over 10 cycles for each c-rate
% electrode = ['n', 'p'];
% for e = electrode
%     for i = 1:length(c_rates)
%         % negative electrode 
%         figure('position', [0 0 600 300]); hold on;  grid on; a = [];
%         results_table = results{i};
%         for j = 1:length(nr_dims)
%             col_name = 'cs' + string(e) + '_nr' + string(nr_dims(j));
%             x = 1:(t_duration(i)*2*10);    
%             y = results_table.(col_name);
%             a = [a, plot(x, y, 'linewidth', 1)];
%         end
%         col_name = 'cs' + string(e) + '_nr5';
%         yline(max(results_table.(col_name)(1:(2*t_duration(i)))), 'r');
%         yline(min(results_table.(col_name)(1:(2*t_duration(i)))), 'r');
%         M1 = "Nr = 5"; M2 = "Nr = 10"; M3 = "Nr = 100";
%         legend(a, [M1, M2, M3]);
%         xlabel('time step (sec)','Fontsize',12,'interpreter','latex')
%         ylabel('cs'+string(e)+' voltage @ c-rate of '+string(c_rates(i)),'Fontsize',12,'interpreter','latex')
%         set(gca,'Fontsize',12);
%         filename = strcat('voltage_over_cycles', '_cs', string(e), '_cr', string(c_rates(i)), '_plot.png');
%         saveas(gcf,fullfile('output', filename)) 
%     end
% end
