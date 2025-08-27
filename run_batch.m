%% 1. download data

% from OSF: https://doi.org/10.17605/OSF.IO/5F6NB

%% 2. compute HFB

% run scripts in order
hfb_analysis;
hfb_ga;

%% 3. run models on HFB

% add path to subfunctions
addpath(fullfile(pwd, 'subfunctions'));

hfb_stat_lme_trl;

% reliability models
for num_trl = 1:20
    hfb_stat_lme_trl_reliability(num_trl);
    for num_sbj = 1:18
        hfb_stat_lme_trl_reliability_sbj(num_sbj, num_trl);
    end
    hfb_stat_lme_trl_reliability_ind_ch(num_trl);
end
        