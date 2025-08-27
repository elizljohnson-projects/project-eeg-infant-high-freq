function hfb_stat_lme_trl_reliability_ind_ch(num_trl)
% HFB_STAT_LME_TRL_RELIABILITY_IND_CH - run linear mixed-effects model on 
% per-subject HFB sleep vs. wake data, using the channel-wise mean based on 
% group-level results, varying the number of trials and repeating the 
% analysis 10 times.
%
% Inputs:
% num_trl = number of trials (e.g., 5)
%
% Example:
% hfb_stat_lme_trl_reliability_ind_ch(5)
%
% Copyright (c) 2023
% EL Johnson, PhD & AM Holubecki, PhD candidate

% set directories
pth = pwd;
datdir = fullfile(pth, 'All_Data_HFB'); % ouput of hfb_ga
savdir = fullfile(datdir, 'reliability', 'individual');
mkdir(savdir);

num_x = 10; % # replications

% load data
load(fullfile(datdir, 'hfb_ga_trl'), 'sleep', 'awake');

% set channels
[ch, ch_idx] = intersect(sleep.label, {'Fz', 'Cz', 'Pz', 'C3', 'C4', 'O1', 'O2'});

% loop through replications
x_rand_s = cell(length(sleep.powspctrm), num_x);
x_rand_a = x_rand_s;
for x = 1:num_x
    for s = 1:length(sleep.powspctrm)
        % randomly sample trials without replacement
        x_rand_s{s,x} = randsample(size(sleep.powspctrm{s},1), num_trl, false);
        x_rand_a{s,x} = randsample(size(awake.powspctrm{s},1), num_trl, false);
    end
    
    % loop through subjects, running the model per subject
    for s = 1:length(sleep.powspctrm)
        % initialize structure for model outputs in FieldTrip-like format
        lme = [];
        lme.dimord = 'chan_freq';
        lme.label = ch';
        lme.freq = sleep.freq;
        
        % extract data for table
        dat = [nanmean(sleep.powspctrm{s}(x_rand_s{s,x},ch_idx),2); ...
            nanmean(awake.powspctrm{s}(x_rand_a{s,x},ch_idx),2)];
        sleep_wake = [zeros(size(sleep.powspctrm{s}(x_rand_s{s,x},1))); ...
            ones(size(awake.powspctrm{s}(x_rand_a{s,x},1)))];
        trl = 1:length(dat);
        trl = trl';
        
        % put in table
        data = table(num2str(trl), num2str(sleep_wake), dat, 'VariableNames', ...
            {'trl', 'awake_sleep', 'dat'});
        
        % run model
        tmp = fitlme(data, 'dat ~ awake_sleep + (1|trl)');
        
        % extract model outputs
        lme.t = tmp.Coefficients{2,4};
        lme.p = tmp.Coefficients{2,6};
        lme.df = tmp.Coefficients{2,5};
        
        % save
        save(fullfile(savdir, ['hfb_stat_lme_s' num2str(s) '_' num2str(num_trl) 'trl_r' num2str(x)]), 'lme');

        clear data dat sleep_wake trl lme
    end
end

end
