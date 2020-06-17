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

%% read the data list: 3T
fdata = [proj_dir, '/data/volbrain_replicability_7T_3R-BRAIN.csv'];
table_vol = readtable(fdata);
%demographic info
subid = table_vol.subid;
age = table_vol.age;
sessid = table_vol.sesid;
scanner = table_vol.scanner;
sex = table_vol.gender; % 1 is male
icv = table_vol.TissueICCm3;
%subcortical volume
Thalamus_L = table_vol.ThalamusLeftCm3;
Thalamus_R = table_vol.ThalamusRightCm3;
Accumbens_L = table_vol.AccumbensLeftCm3;
Accumbens_R = table_vol.AccumbensRightCm3;
Amygdala_L = table_vol.AmygdalaLeftCm3;
Amygdala_R = table_vol.AmygdalaRightCm3;
Caudate_L = table_vol.CaudateLeftCm3;
Caudate_R = table_vol.CaudateRightCm3;
Pallidum_L = table_vol.GlobusPallidusLeftCm3;
Pallidum_R = table_vol.GlobusPallidusRightCm3;
Putamen_L = table_vol.PutamenLeftCm3;
Putamen_R = table_vol.PutamenRightCm3;
Hippocampus_L = table_vol.HippocampusLeftCm3;
Hippocampus_R = table_vol.HippocampusRightCm3;

%% Estimate Laterality for Subcortical Regions
names_subcort = {'Thalamus', 'Caudate', 'Putamen', 'Pallidum', ...
    'Hippocampus', 'Amygdala', 'Accumbens'};
%estimate asymmetry index (AI)
vol_u103_subcort_L = [Thalamus_L Caudate_L Putamen_L Pallidum_L ...
    Hippocampus_L Amygdala_L Accumbens_L];
vol_u103_subcort_R = [Thalamus_R Caudate_R Putamen_R Pallidum_R ...
    Hippocampus_R Amygdala_R Accumbens_R];
LvR_subcort = vol_u103_subcort_L - vol_u103_subcort_R;
LaR_subcort = vol_u103_subcort_L + vol_u103_subcort_R;
AI_u103_subcort = 100*(2*LvR_subcort./LaR_subcort);

%% Paired stats tests on laterality: GE750 (just for test)
idx_GE = contains(scanner,'IPCAS_MR750');
hGE = zeros(7,1); pGE = ones(7,1);
ciGE = zeros(7,2); statsGE = cell(7,1);
for idx=1:7
    tmpL = vol_u103_subcort_L(idx_GE,idx);
    tmpR = vol_u103_subcort_R(idx_GE,idx);
    [hGE(idx),pGE(idx),ciGE(idx,:),statsGE{idx,1}] = ...
        ttest(tmpL,tmpR,'Alpha',0.05/7);
end

%% Paired stats tests on laterality: SMS just for test
idx_SMS = contains(scanner,'IBP_Prisma');
hSMS = zeros(7,1); pSMS = ones(7,1);
ciSMS = zeros(7,2); statsSMS = cell(7,1);
for idx=1:7
    tmpL = vol_u103_subcort_L(idx_SMS,idx);
    tmpR = vol_u103_subcort_R(idx_SMS,idx);
    [hSMS(idx),pSMS(idx),ciSMS(idx,:),statsSMS{idx,1}] = ...
        ttest(tmpL,tmpR,'Alpha',0.05/7);
end

%% Paired stats tests on laterality: SMS7T
idx_SMS7T = contains(scanner,'IBP_Magnetom');
hSMS7T = zeros(7,1); pSMS7T = ones(7,1);
ciSMS7T = zeros(7,2); statsSMS7T = cell(7,1);
effSMS7T = zeros(7,3); pwrSMS7T = zeros(7,1);
for idx=1:7
    tmpL = vol_u103_subcort_L(idx_SMS7T,idx);
    tmpR = vol_u103_subcort_R(idx_SMS7T,idx);
    [hSMS7T(idx),pSMS7T(idx),ciSMS7T(idx,:),tmpstats] = ...
        ccs_core_ttests(tmpL,tmpR,'paired',0.05/7);
    %effect size
    effSMS7T(idx,:) = tmpstats.eftsize;
    pwrSMS7T(idx) = tmpstats.power;
    statsSMS7T{idx,1} = tmpstats;
end
%for plot in prism
[~,idx_p] = sort(pSMS7T); 
effSMS7T_prism = -effSMS7T(idx_p,[3 2 1])';
pwrSMS7T_prism = pwrSMS7T(idx_p);
%save the results
names_subcort_sort = names_subcort(idx_p);
disp(names_subcort(idx_p));
fresults = [proj_dir '/zxn/results/replicability/volbrain_ibp7t_LR_diff.mat'];
save(fresults,'hSMS7T','pSMS7T','ciSMS7T','statsSMS7T','effSMS7T','pwrSMS7T',...
    'names_subcort')

%% PreStat tests effects of age and icv on AI
[rGE_age,pGE_age] = corr(AI_u103_subcort(idx_GE,:),age(idx_GE));
[rGE_icv,pGE_icv] = corr(AI_u103_subcort(idx_GE,:),icv(idx_GE));
[rSMS_age,pSMS_age] = corr(AI_u103_subcort(idx_SMS,:),age(idx_SMS));
[rSMS_icv,pSMS_icv] = corr(AI_u103_subcort(idx_SMS,:),icv(idx_SMS));
[rSMS7T_age,pSMS7T_age] = corr(AI_u103_subcort(idx_SMS7T,:),age(idx_SMS7T));
[rSMS7T_icv,pSMS7T_icv] = corr(AI_u103_subcort(idx_SMS7T,:),icv(idx_SMS7T));

%% Sex-related differences in AI: GE (just for test)
cov_dm = IPN_demean([age(idx_GE) icv(idx_GE)]);
hAI_GE = zeros(7,1); pAI_GE = ones(7,1);
ciAI_GE = zeros(7,2); statsAI_GE = cell(7,1);
tmpsex = sex(idx_GE);
for idx=1:7
    tmpAI = AI_u103_subcort(idx_GE,idx);
    [~,~,tmpAI_adj] = regress(tmpAI,cov_dm);
    %female
    tmpAI_female_adj = tmpAI_adj(tmpsex==0);
    %male
    tmpAI_male_adj = tmpAI_adj(tmpsex==1);
    %two-sample t-tests
    [hAI_GE(idx),pAI_GE(idx),ciAI_GE(idx,:),statsAI_GE{idx,1}] = ...
        ttest2(tmpAI_female_adj,tmpAI_male_adj,'Alpha',0.05/7);
end

%% Sex-related differences in AI: SMS (just for test)
cov_dm = IPN_demean([age(idx_SMS) icv(idx_SMS)]);
hAI_SMS = zeros(7,1); pAI_SMS = ones(7,1);
ciAI_SMS = zeros(7,2); statsAI_SMS = cell(7,1);
tmpsex = sex(idx_SMS);
for idx=1:7
    tmpAI = AI_u103_subcort(idx_SMS,idx);
    [~,~,tmpAI_adj] = regress(tmpAI,cov_dm);
    %female
    tmpAI_female_adj = tmpAI_adj(tmpsex==0);
    %male
    tmpAI_male_adj = tmpAI_adj(tmpsex==1);
    %two-sample t-tests
    [hAI_SMS(idx),pAI_SMS(idx),ciAI_SMS(idx,:),statsAI_SMS{idx,1}] = ...
        ttest2(tmpAI_female_adj,tmpAI_male_adj,'Alpha',0.05/7);
end

%% Sex-related differences in AI: SMS7T
cov_dm = IPN_demean([age(idx_SMS7T) icv(idx_SMS7T)]);
hAI_SMS7T = zeros(7,1); pAI_SMS7T = ones(7,1);
ciAI_SMS7T = zeros(7,2); statsAI_SMS7T = cell(7,1);
effAI_SMS7T = zeros(7,3); pwrAI_SMS7T = zeros(7,1);
tmpsex = sex(idx_SMS7T);
for idx=1:7
    tmpAI = AI_u103_subcort(idx_SMS7T,idx);
    [~,~,tmpAI_adj] = regress(tmpAI,cov_dm);
    %female
    tmpAI_female_adj = tmpAI_adj(tmpsex==0);
    %male
    tmpAI_male_adj = tmpAI_adj(tmpsex==1);
    %two-sample t-tests
    [hAI_SMS7T(idx),pAI_SMS7T(idx),ciAI_SMS7T(idx,:),tmpstats] = ...
        ccs_core_ttests(tmpAI_female_adj,tmpAI_male_adj,'independent',0.05/7);
    %effect size and power
    effAI_SMS7T(idx,:) = tmpstats.eftsize;
    pwrAI_SMS7T(idx) = tmpstats.power;
    statsAI_SMS7T{idx,1} = tmpstats;
end
%for plot in prism
[~,idx_p] = sort(pAI_SMS7T); 
effAI_SMS7T_prism = -effAI_SMS7T(idx_p,[3 2 1])';
pwrAI_SMS7T_prism = pwrAI_SMS7T(idx_p);
%save the results
disp(names_subcort(idx_p));
fresults = [proj_dir '/zxn/results/replicability/volbrain_ibp7t_sexAI_diff.mat'];
save(fresults,'hAI_SMS7T','pAI_SMS7T','ciAI_SMS7T','statsAI_SMS7T',...
    'effAI_SMS7T','pwrAI_SMS7T','names_subcort')
