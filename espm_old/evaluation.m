%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Evaluate stability & accuracy of model %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
warning off
clc

scheme = 'fdm';
ref_scheme = 'espm_fdm';

%% Reference solutions
% load data from previous runs
ref_iapp24 = load(fullfile('output',append(ref_scheme,'_nr1000_cr1_cyc0.mat')));
ref_iapp48 = load(fullfile('output',append(ref_scheme,'_nr1000_cr2_cyc0.mat')));
ref_iapp96 = load(fullfile('output',append(ref_scheme,'_nr1000_cr3_cyc0.mat')));
ref_data = [ref_iapp24 ref_iapp48 ref_iapp96];

ref_table = table;
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

    ref_table.i = [
        ref.all_data.V(end); 
        ref.all_data.csn(end,end);
        ref.all_data.csp(end,end);
        dot(weights_n, ref.all_data.csn(:,end));
        dot(weights_p, ref.all_data.csp(:,end));
    ];
end

ref_table.Properties.RowNames = {'voltage';'cathode_surf_con';'anode_surf_con';'cathode_avg_con';'anode_avg_con'};
ref_table.Properties.VariableNames = {'iapp24'}; %'iapp48';'iapp96';
writetable(ref_table, fullfile('ref_table.dat'))


%% Experiments to evaluate accuracy
% load data from previous runs
ref_iapp24 = load(append(scheme,'_nr1000_cr1_cyc0.mat'));
ref_iapp48 = load(append(scheme,'_nr1000_cr2_cyc0.mat'));
ref_iapp96 = load(append(scheme,'_nr1000_cr3_cyc0.mat'));
ref_data = [ref_iapp24 ref_iapp48 ref_iapp96];

%% Plot V_RMS, Csurf_RMS, Cavg_RMS