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

%% read the subject list
flist_u100 = [proj_dir, '/data/behavioral_unrelated100_hcp.csv'];
table_u100_BH = readtable(flist_u100);
subjects_id = table_u100_BH.Subject;

%% Load FreeSurfer table - HCP1200 release
flist_1200 = [proj_dir '/data/behavioral_s1200_hcp.csv'];
table_1200_BH = readtable(flist_1200);
subjects_1200 = table_1200_BH.Subject;
%1200 restricted info
flist_1200R = [proj_dir '/data/hcp_restricted.csv']; %get this info from HCP
table_1200R_BH = readtable(flist_1200R);
subjects_1200R = table_1200R_BH.Subject; 
%get demographic info
[~, subID_3T_BH, ~] = intersect(subjects_1200,subjects_id);
sex = table_1200_BH.Gender(subID_3T_BH);
[~, subID_3T_RBH, ~] = intersect(subjects_1200R,subjects_id);
age = table_1200R_BH.Age_in_Yrs(subID_3T_RBH);
hand = table_1200R_BH.Handedness(subID_3T_RBH);
%get freesurfer subcortical volumes
fTable = [proj_dir '/data/freesurfer_hcp.csv'];
table_1200_FS = readtable(fTable);
subjects_FS = table_1200_FS.Subject; 
[~, subID_3T_FS, ~] = intersect(subjects_FS,subjects_id);
table_u100 = table_1200_FS(subID_3T_FS,:);

%% Estimate Laterality for Subcortical Regions
ICV_u100 = table_u100.FS_IntraCranial_Vol;
idx_subcort = [25:28 32 33 35 43:49];
names_subcort = table_u100.Properties.VariableNames(idx_subcort);
table_u100_subcort = table_u100(:,idx_subcort);
vol_u100_subcort_L = table2array(table_u100_subcort(:,1:7));
vol_u100_subcort_R = table2array(table_u100_subcort(:,8:14));
LvR_subcort = vol_u100_subcort_L - vol_u100_subcort_R;
LaR_subcort = vol_u100_subcort_L + vol_u100_subcort_R;
AI_u100_subcort = 100*(2*LvR_subcort./LaR_subcort);

%% Paired stats tests on laterality
h = zeros(7,1); p = ones(7,1);
ci = zeros(7,2); stats = cell(7,1);
eff = zeros(7,3); pwr = zeros(7,1);
for idx=1:7
    tmpL = vol_u100_subcort_L(:,idx);
    tmpR = vol_u100_subcort_R(:,idx);
    %paired t-test
    [h(idx),p(idx),ci(idx,:),tmpstats] = ...
        ccs_core_ttests(tmpL,tmpR,'paired',0.05/7);
    %effect size
    eff(idx,:) = tmpstats.eftsize;
    pwr(idx) = tmpstats.power;
    stats{idx,1} = tmpstats;
end
%for plot in prism
[p_prism,idx_p] = sort(p);
eff_prism = -eff(idx_p,[3 2 1])';
pwr_prism = pwr(idx_p);
%save the results
disp(names_subcort(idx_p));
fresults = [proj_dir '/zxn/results/replicability/volbrain_hcpu100_LR_diff.mat'];
save(fresults,'h','p','ci','stats','eff','pwr','names_subcort')

%% PreStat tests effects of age,sex,hand and icv on AI
[r_age,p_age] = corr(AI_u100_subcort,age);
[r_hand,p_hand] = corr(AI_u100_subcort,hand);
[r_icv,p_icv] = corr(AI_u100_subcort,ICV_u100);

%% Sex-related differences in AI
idx_female = contains(sex,'F');
idx_male = contains(sex,'M');
%regress covariates seperately for males and females
hAI = zeros(7,1); pAI = ones(7,1);
ciAI = zeros(7,2); statsAI = cell(7,1);
effAI = zeros(7,3); pwrAI = zeros(7,1);
%regress covariates
cov_dm = IPN_demean([age hand ICV_u100]);
for idx=1:7
    tmpAI = AI_u100_subcort(:,idx);
    [~,~,tmpAI_adj] = regress(tmpAI,cov_dm);
    %female
    tmpAI_female_adj = tmpAI_adj(idx_female);
    %male
    tmpAI_male_adj = tmpAI_adj(idx_male);
    %two-sample t-tests
    [hAI(idx),pAI(idx),ciAI(idx,:),tmpstats] = ...
        ccs_core_ttests(tmpAI_female_adj,tmpAI_male_adj,'independent',0.05/7);
    %effect size and power
    effAI(idx,:) = tmpstats.eftsize;
    pwrAI(idx) = tmpstats.power;
    statsAI{idx,1} = tmpstats;
end
%for plot in prism
[~,idx_p] = sort(pAI); 
effAI_prism = -effAI(idx_p,[3 2 1])';
pwrAI_prism = pwrAI(idx_p);
%save the results
disp(names_subcort(idx_p));
fresults = [proj_dir '/zxn/results/replicability/freesurfer_hcpu100_sexAI_diff.mat'];
save(fresults,'hAI','pAI','ciAI','statsAI','effAI','pwrAI','names_subcort')

%% re-examine the whole hcp sample
[~, subID_FS_BH, ~] = intersect(subjects_1200,subjects_FS);
sex_FS = table_1200_BH.Gender(subID_FS_BH);
[~, subID_FS_RBH, ~] = intersect(subjects_1200R,subjects_FS);
age_FS = table_1200R_BH.Age_in_Yrs(subID_FS_RBH);
hand_FS = table_1200R_BH.Handedness(subID_FS_RBH);
% estimate Laterality for Subcortical Regions
ICV = table_1200_FS.FS_IntraCranial_Vol;
table_1200_subcort = table_1200_FS(:,idx_subcort);
vol_1200_subcort_L = table2array(table_1200_subcort(:,1:7));
vol_1200_subcort_R = table2array(table_1200_subcort(:,8:14));
LvR_subcort = vol_1200_subcort_L - vol_1200_subcort_R;
LaR_subcort = vol_1200_subcort_L + vol_1200_subcort_R;
AI_1200_subcort = 100*(2*LvR_subcort./LaR_subcort);
% paired stats tests on laterality
hFS = zeros(7,1); pFS = ones(7,1);
ciFS = zeros(7,2); statsFS = cell(7,1);
effFS = zeros(7,3); pwrFS = zeros(7,1);
for idx=1:7
    tmpL = vol_1200_subcort_L(:,idx);
    tmpR = vol_1200_subcort_R(:,idx);
    [hFS(idx),pFS(idx),ciFS(idx,:),tmpstats] = ...
        ccs_core_ttests(tmpL,tmpR,'paired',0.05/7);
    %effect size
    effFS(idx,:) = tmpstats.eftsize;
    pwrFS(idx) = tmpstats.power;
    statsFS{idx,1} = tmpstats;
end
%for plot in prism
[pFS_prism,idx_p] = sort(pFS);
effFS_prism = -effFS(idx_p,[3 2 1])';
pwrFS_prism = pwrFS(idx_p);
%save the results
disp(names_subcort(idx_p));
fresults = [proj_dir '/zxn/results/replicability/volbrain_hcp1113_LR_diff.mat'];
save(fresults,'hFS','pFS','ciFS','statsFS','effFS','pwrFS','names_subcort')
% preStat tests effects of age,sex,hand and icv on AI
[rFS_age,pFS_age] = corr(AI_1200_subcort,age_FS);
[rFS_hand,pFS_hand] = corr(AI_1200_subcort,hand_FS);
[rFS_icv,pFS_icv] = corr(AI_1200_subcort,ICV);
% sex-related differences in AI
idx_female = contains(sex_FS,'F');
idx_male = contains(sex_FS,'M');
%regress covariates seperately for males and females
hAI_FS = zeros(7,1); pAI_FS = ones(7,1);
ciAI_FS = zeros(7,2); statsAI_FS = cell(7,1);
effAI_FS = zeros(7,3); pwrAI_FS = zeros(7,1);
%regress covariates
cov_dm = IPN_demean([age_FS hand_FS ICV]);
for idx=1:7
    tmpAI = AI_1200_subcort(:,idx);
    [~,~,tmpAI_adj] = regress(tmpAI,cov_dm);
    %female
    tmpAI_female_adj = tmpAI_adj(idx_female);
    %male
    tmpAI_male_adj = tmpAI_adj(idx_male);
    %two-sample t-tests
    [hAI_FS(idx),pAI_FS(idx),ciAI_FS(idx,:),tmpstats] = ...
        ccs_core_ttests(tmpAI_female_adj,tmpAI_male_adj,'independent',0.05/7);
    %effect size and power
    effAI_FS(idx,:) = tmpstats.eftsize;
    pwrAI_FS(idx) = tmpstats.power;
    statsAI_FS{idx,1} = tmpstats;
end
%for plot in prism
[~,idx_p] = sort(pAI_FS); 
effAI_FS_prism = -effAI_FS(idx_p,[3 2 1])';
pwrAI_FS_prism = pwrAI_FS(idx_p);
%save the results
disp(names_subcort(idx_p));
fresults = [proj_dir '/zxn/results/replicability/freesurfer_hcpu1113_sexAI_diff.mat'];
save(fresults,'hAI_FS','pAI_FS','ciAI_FS','statsAI_FS',...
    'effAI_FS','pwrAI_FS','names_subcort')
