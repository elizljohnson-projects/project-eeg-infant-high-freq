% HFB_STAT_LME_TRL - run linear mixed-effects model on HFB sleep vs. wake
% data, including all subjects and all trials per subject. Plot results on 
% the scalp topography. Outputs are FDR-corrected across channels using the 
% subfunction fdr.m.
%
% Ensure FieldTrip is correcty added to the MATLAB path:
%   addpath <path to fieldtrip home directory>
%   ft_defaults
%
% Copyright (c) 2023
% EL Johnson, PhD & AM Holubecki, PhD candidate

% set directories
pth = pwd;
datdir = fullfile(pth, 'All_Data_HFB'); % ouput of hfb_ga
savdir = datdir;

% load data
load(fullfile(datdir, 'hfb_ga_trl'), 'sleep', 'awake');

% initialize structure for model outputs in FieldTrip-like format
lme = [];
lme.dimord = 'chan_freq';
lme.label = sleep.label;
lme.freq = sleep.freq;
lme.t = nan(length(lme.label), 1);
lme.p = lme.t;
lme.df = lme.t;

% loop through channels, running the model per channel
hold = [];
for e = 1:length(sleep.label)
    % extract data for table
    sid = [];
    trl = [];
    sleep_wake = [];
    dat = [];
    for s = 1:length(sleep.powspctrm)
        tmp_dat = [sleep.powspctrm{s}(:,e); awake.powspctrm{s}(:,e)];
        tmp_s_a = [zeros(size(sleep.powspctrm{s}(:,e))); ...
            ones(size(awake.powspctrm{s}(:,e)))];
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
    hold = cat(2, hold, dat);
    
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

% put all data in master table for later use
data = table(sid, trl, sleep_wake, hold(:,1), hold(:,2), hold(:,3), hold(:,4), ...
    hold(:,5), hold(:,6), hold(:,7), hold(:,8), hold(:,9), hold(:,10), hold(:,11), ...
    hold(:,12), hold(:,13), hold(:,14), hold(:,15), hold(:,16), hold(:,17), ...
    hold(:,18), hold(:,19), 'VariableNames', {'sid', 'trl', 'sleep_wake', lme.label{1}, ...
    lme.label{2}, lme.label{3}, lme.label{4}, lme.label{5}, lme.label{6}, ...
    lme.label{7}, lme.label{8}, lme.label{9}, lme.label{10}, lme.label{11}, ...
    lme.label{12}, lme.label{13}, lme.label{14}, lme.label{15}, lme.label{16}, ...
    lme.label{17}, lme.label{18}, lme.label{19}});

% save
save(fullfile(savdir, 'hfb_stat_lme_trl'), 'lme', 'data');

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

% save figure
print(fullfile(savdir, 'hfb_stat_lme_trl'), '-dtiff');
