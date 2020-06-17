%% dir settings (may not usable for you and you have to change them...)
clear all; clc
tmplist = dir; tmpdir = tmplist.folder;
proj_dir = tmpdir(1:end-9);
ccs_dir = '/Users/mac/Projects/CCS'; %you may need to change this
ccs_matlab = [ccs_dir '/matlab'];
ccs_vistool = [ccs_dir '/vistool'];

%% Set up the path to matlab function in Freesurfer release
addpath(genpath(ccs_matlab)) %ccs matlab scripts
addpath(genpath(ccs_vistool)) %ccs matlab scripts

%% read the results list: 3T GE750
fdata = [proj_dir, '/wys/results/reproducibility/', ...
    'volbrain_ipcas_for-3T_LR_diff_by_wys.csv'];
table_ptest = readtable(fdata);
table_ptest_subcort = table_ptest(:,[1:4 11 12 14 15]);
p_subcort = table2array(table_ptest_subcort(2,2:end));
[~,idx_p] = sort(p_subcort); p_subcort(idx_p)*7 %FDR corrected
names_subcort = table_ptest_subcort(:,2:end).Properties.VariableNames;
names_subcort = names_subcort(idx_p);
stats_subcort = table2array(table_ptest_subcort(:,idx_p+1));
effsize_subcort = stats_subcort([7 6 8],:);
effsize_subcort_flip = -effsize_subcort;

%% read the results list: matching to subjects with 7T SMS
fdata = [proj_dir, '/wys/results/reproducibility/', ...
    'volbrain_ipcas_for-7T_LR_diff_by_wys.csv'];
table_ptest_7T = readtable(fdata);
table_ptest_subcort_7T = table_ptest_7T(:,[1:4 11 12 14 15]);
p_subcort_7T = table2array(table_ptest_subcort_7T(2,2:end));
[~,idx_p] = sort(p_subcort_7T); p_subcort_7T(idx_p)*7 %FDR corrected
names_subcort_7T = table_ptest_subcort_7T(:,2:end).Properties.VariableNames;
names_subcort_7T = names_subcort_7T(idx_p);
stats_subcort_7T = table2array(table_ptest_subcort_7T(:,idx_p+1));
effsize_subcort_7T = stats_subcort_7T([7 6 8],:);
effsize_subcort_7T_flip = - effsize_subcort_7T;

%% read the results list: 3T SMS
fdata = [proj_dir, '/wys/results/replicability/', ...
    'volbrain_ibp3t_LR_diff_by_wys.csv'];
table_ptest_ibp3t = readtable(fdata);
table_ptest_subcort_ibp3t = table_ptest_ibp3t(:,[1:4 11 12 14 15]);
p_subcort_ibp3t = table2array(table_ptest_subcort_ibp3t(2,2:end));
[~,idx_p_ibp3t] = sort(p_subcort_ibp3t); p_subcort_ibp3t(idx_p_ibp3t)*7 %FDR corrected
names_subcort_ibp3t = table_ptest_subcort_ibp3t(:,2:end).Properties.VariableNames;
names_subcort_ibp3t = names_subcort_ibp3t(idx_p_ibp3t);
stats_subcort_ibp3t = table2array(table_ptest_subcort_ibp3t(:,idx_p_ibp3t+1));
effsize_subcort_ibp3t = stats_subcort_ibp3t([7 6 8],:);
effsize_subcort_ibp3t_flip = -effsize_subcort_ibp3t;

%% read the results list: 7T SMS - real
fdata = [proj_dir, '/wys/results/replicability/', ...
    'volbrain_ibp7t_LR_diff_by_wys.csv'];
table_ptest_ibp7t = readtable(fdata);
table_ptest_subcort_ibp7t = table_ptest_ibp7t(:,[1:4 11 12 14 15]);
p_subcort_ibp7t = table2array(table_ptest_subcort_ibp7t(2,2:end));
[~,idx_p_ibp7t] = sort(p_subcort_ibp7t); p_subcort_ibp7t(idx_p_ibp7t)*7 %FDR corrected
names_subcort_ibp7t = table_ptest_subcort_ibp7t(:,2:end).Properties.VariableNames;
names_subcort_ibp7t = names_subcort_ibp7t(idx_p_ibp7t);
stats_subcort_ibp7t = table2array(table_ptest_subcort_ibp7t(:,idx_p_ibp7t+1));
effsize_subcort_ibp7t = stats_subcort_ibp7t([7 6 8],:);
effsize_subcort_ibp7t_flip = -effsize_subcort_ibp7t;

%% read the results list: ICC
fdata = [proj_dir, '/wys/results/reliability/volbrain_icc_updated.csv'];
table_reliability = readtable(fdata);
%index for subcortical regions
idx_subcort = [3 7 12 36 41 49 56];
%extract variabilities
table_var_subcortL = table_reliability(6:12,idx_subcort);
table_var_subcortR = table_reliability(6:12,idx_subcort+1);
table_var_subcortAI = table_reliability(6:12,idx_subcort-1);
%extract center icc
table_iccC_subcortL = table_reliability(1:3,idx_subcort);
table_iccC_subcortR = table_reliability(1:3,idx_subcort+1);
table_iccC_subcortAI = table_reliability(1:3,idx_subcort-1);
%extract low-bound icc
fdata = [proj_dir, '/wys/results/reliability/volbrain_icc_ci25.csv'];
table_iccL = readtable(fdata);
table_iccL_subcortL = table_iccL(1:3,idx_subcort);
table_iccL_subcortR = table_iccL(1:3,idx_subcort+1);
table_iccL_subcortAI = table_iccL(1:3,idx_subcort-1);
%extract up-bound icc
fdata = [proj_dir, '/wys/results/reliability/volbrain_icc_ci975.csv'];
table_iccU = readtable(fdata);
table_iccU_subcortL = table_iccU(1:3,idx_subcort);
table_iccU_subcortR = table_iccU(1:3,idx_subcort+1);
table_iccU_subcortAI = table_iccU(1:3,idx_subcort-1);

%% prepare for plots in prism
names_subcort = table_reliability(:,idx_subcort+2).Properties.VariableNames;
grpm = table2array(table_reliability(4,idx_subcort+2));
[~,idx_subcort_update] = sort(grpm,'descend');
disp(names_subcort(idx_subcort_update))
%variabilities: L and R
table_varSubj_subcort = [-table2array(table_var_subcortL(1,idx_subcort_update))' ...
    table2array(table_var_subcortR(1,idx_subcort_update))'];
table_varSite_subcort = [-(table2array(table_var_subcortL(2,idx_subcort_update))' + ...
    table2array(table_var_subcortL(5,idx_subcort_update))' + ...
    table2array(table_var_subcortL(7,idx_subcort_update))') ...
    (table2array(table_var_subcortR(2,idx_subcort_update))' + ...
    table2array(table_var_subcortR(5,idx_subcort_update))' + ...
    table2array(table_var_subcortR(7,idx_subcort_update))')];
table_varSess_subcort = [-(table2array(table_var_subcortL(3,idx_subcort_update))' + ...
    table2array(table_var_subcortL(6,idx_subcort_update))' + ...
    table2array(table_var_subcortL(7,idx_subcort_update))') ...
    (table2array(table_var_subcortR(3,idx_subcort_update))' + ...
    table2array(table_var_subcortR(6,idx_subcort_update))' + ...
    table2array(table_var_subcortR(7,idx_subcort_update))')]; 
%variabilities: L vs R
table_varSubj_subcortAI = table2array(table_var_subcortAI(1,idx_subcort_update))';
table_varSite_subcortAI = table2array(table_var_subcortAI(2,idx_subcort_update))' + ...
    table2array(table_var_subcortL(5,idx_subcort_update))' + ...
    table2array(table_var_subcortL(7,idx_subcort_update))';
table_varSess_subcortAI = table2array(table_var_subcortL(3,idx_subcort_update))' + ...
    table2array(table_var_subcortL(6,idx_subcort_update))' + ...
    table2array(table_var_subcortL(7,idx_subcort_update))';
%overall icc
table_icc_subcortL = [table2array(table_iccU_subcortL(3,idx_subcort_update)); ...
    table2array(table_iccC_subcortL(3,idx_subcort_update)); ...
    table2array(table_iccL_subcortL(3,idx_subcort_update))];
table_icc_subcortR = [table2array(table_iccU_subcortR(3,idx_subcort_update)); ...
    table2array(table_iccC_subcortR(3,idx_subcort_update)); ...
    table2array(table_iccL_subcortR(3,idx_subcort_update))];
table_icc_subcortAI = [table2array(table_iccU_subcortAI(3,idx_subcort_update)); ...
    table2array(table_iccC_subcortAI(3,idx_subcort_update)); ...
    table2array(table_iccL_subcortAI(3,idx_subcort_update))];
%test-retest icc
table_iccSess_subcortL = [table2array(table_iccU_subcortL(1,idx_subcort_update)); ...
    table2array(table_iccC_subcortL(1,idx_subcort_update)); ...
    table2array(table_iccL_subcortL(1,idx_subcort_update))];
table_iccSess_subcortR = [table2array(table_iccU_subcortR(1,idx_subcort_update)); ...
    table2array(table_iccC_subcortR(1,idx_subcort_update)); ...
    table2array(table_iccL_subcortR(1,idx_subcort_update))];
table_iccSess_subcortAI = [table2array(table_iccU_subcortAI(1,idx_subcort_update)); ...
    table2array(table_iccC_subcortAI(1,idx_subcort_update)); ...
    table2array(table_iccL_subcortAI(1,idx_subcort_update))];
%inter-site icc
table_iccSite_subcortL = [table2array(table_iccU_subcortL(2,idx_subcort_update)); ...
    table2array(table_iccC_subcortL(2,idx_subcort_update)); ...
    table2array(table_iccL_subcortL(2,idx_subcort_update))];
table_iccSite_subcortR = [table2array(table_iccU_subcortR(2,idx_subcort_update)); ...
    table2array(table_iccC_subcortR(2,idx_subcort_update)); ...
    table2array(table_iccL_subcortR(2,idx_subcort_update))];
table_iccSite_subcortAI = [table2array(table_iccU_subcortAI(2,idx_subcort_update)); ...
    table2array(table_iccC_subcortAI(2,idx_subcort_update)); ...
    table2array(table_iccL_subcortAI(2,idx_subcort_update))];