function hfb_stat_lme_trl_reliability(num_trl)
% HFB_STAT_LME_TRL_RELIABILITY - run linear mixed-effects model on HFB 
% sleep vs. wake data, including all subjects, varying the number of trials
% per subject, and repeating the analysis 10 times. Plot results on the 
% scalp topography. Outputs are FDR-corrected across channels using the 
% subfunction fdr.m.
%
% Ensure FieldTrip is correcty added to the MATLAB path:
%   addpath <path to fieldtrip home directory>
%   ft_defaults
%
% Inputs:
% num_trl = number of trials (e.g., 5)
%
% Example:
% hfb_stat_lme_trl_reliability(5)
%
% Copyright (c) 2023
% EL Johnson, PhD & AM Holubecki, PhD candidate

% set directories
pth = pwd;
datdir = fullfile(pth, 'All_Data_HFB'); % ouput of hfb_ga
savdir = fullfile(datdir, 'reliability');
mkdir(savdir);

num_x = 10; % # replications

% load data
load(fullfile(datdir, 'hfb_ga_trl'), 'sleep', 'awake');

% loop through replications
x_rand_s = cell(length(sleep.powspctrm), num_x);
x_rand_a = x_rand_s;
for x = 1:num_x
    % randomly sample trials without replacement
    for s = 1:length(sleep.powspctrm)
        x_rand_s{s,x} = randsample(size(sleep.powspctrm{s},1), num_trl, false);
        x_rand_a{s,x} = randsample(size(awake.powspctrm{s},1), num_trl, false);
    end
    
    % initialize structure for model outputs in FieldTrip-like format
    lme = [];
    lme.dimord = 'chan_freq';
    lme.label = sleep.label;
    lme.freq = sleep.freq;
    lme.t = nan(length(lme.label), 1);
    lme.p = lme.t;
    lme.df = lme.t;
    
    % loop through channels, running the model per channel
    for e = 1:length(sleep.label)
        % extract data for table
        sid = [];
        trl = [];
        sleep_wake = [];
        dat = [];
        for s = 1:length(sleep.powspctrm)
            tmp_dat = [sleep.powspctrm{s}(x_rand_s{s,x},e); ...
                awake.powspctrm{s}(x_rand_a{s,x},e)];
            tmp_s_a = [zeros(size(sleep.powspctrm{s}(x_rand_s{s,x},e))); ...
                ones(size(awake.powspctrm{s}(x_rand_a{s,x},e)))];
            tmp_sid = repmat(s, size(tmp_dat));
            tmp_trl = 1:length(tmp_dat);
            tmp_trl = tmp_trl';
    
            dat = cat(1, dat, tmp_dat);
            sleep_wake = cat(1, sleep_wake, tmp_s_a);
            sid = cat(1, sid, tmp_sid);
            trl = cat(1, trl, tmp_trl);
            clear tmp*
        end

        % put in table
        data = table(num2str(sid), num2str(trl), num2str(sleep_wake), dat, ...
            'VariableNames', {'sid', 'trl', 'awake_sleep', 'dat'});
        
        % run model
        tmp = fitlme(data, 'dat ~ awake_sleep + (1|sid) + (1|sid:trl)');
    
        % extract model outputs
        lme.t(e) = tmp.Coefficients{2,4};
        lme.p(e) = tmp.Coefficients{2,6};
        lme.df(e) = tmp.Coefficients{2,5};

        clear data tmp
    end

    % apply FDR correction
    [~, lme.fdr_th] = fdr(lme.p, 0.05);
    lme.fdr = lme.p < lme.fdr_th;
    
    % save
    save(fullfile(savdir, ['hfb_stat_lme_' num2str(num_trl) 'trl_r' num2str(x)]), 'lme');

    % change clinical to experimental labels so FieldTrip recognizes them
    [~,idx] = intersect(lme.label, 'T3'); lme.label{idx} = 'T7';
    [~,idx] = intersect(lme.label, 'T4'); lme.label{idx} = 'T8';
    [~,idx] = intersect(lme.label, 'T5'); lme.label{idx} = 'P7';
    [~,idx] = intersect(lme.label, 'T6'); lme.label{idx} = 'P8';
    
    % plot results
    cfg = [];
    cfg.layout = 'elec1020.lay';
    cfg.colorbar = 'yes';
    cfg.gridscale = 300; % high resolution
    cfg.shading = 'interp';
    cfg.style = 'straight';
    cfg.highlight = 'on';
    cfg.highlightsymbol = '.';
    cfg.highlightsize = 30;
    cfg.parameter = 't';
    cfg.zlim = 'maxabs';
    cfg.colormap = '-RdGy';
    cfg.highlightchannel = lme.label(lme.fdr); % mark significant channels
    
    figure; ft_topoplotER(cfg, lme);

    clear lme
    
    % save figure
    print(fullfile(savdir, ['hfb_stat_lme_' num2str(num_trl) 'trl_r' num2str(x)]), '-dtiff');
end

end
