% HFB_GA - compile HFB data across subjects, maintaining the trial 
% dimension.
%
% Copyright (c) 2023
% EL Johnson, PhD & AM Holubecki, PhD candidate

% set directories
pth = pwd;
datdir = fullfile(pth, 'All_Data_HFB'); % ouput of hfb_analysis
savdir = datdir;

% list of channels in the standard clinical EEG headcap
ch_all = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2','F7','F8',...
    'T3','T4','T5','T6','Fz','Cz','Pz'};

% get list of datasets
files = dir(fullfile(datdir, 'Subj*'));

% initialize FieldTrip-like grand average structures
sleep = [];
sleep.dimord = 'subj_rpt_chan_freq';
sleep.label = ch_all;
sleep.powspctrm = cell(length(files), 1);
awake = sleep;

% loop through datasets
for s = 1:length(files)
    % load data
    load(fullfile(datdir, files(s).name));
    
    % finish initializing structures
    if s == 1
        sleep.freq = sleep_eeg.freq;
        awake.freq = sleep.freq;
    end
    
    % populate structures
    sleep.powspctrm{s} = nan(size(sleep_eeg.powspctrm,1), length(ch_all));
    [~,ch_sleep] = intersect(ch_all, sleep_eeg.label, 'stable');
    sleep.powspctrm{s}(:,ch_sleep) = sleep_eeg.powspctrm;

    awake.powspctrm{s} = nan(size(awake_eeg.powspctrm,1), length(ch_all));
    [~,ch_awake] = intersect(ch_all, awake_eeg.label, 'stable');
    awake.powspctrm{s}(:,ch_awake) = awake_eeg.powspctrm;
    
    clear *eeg ch_sleep ch_awake
end

% save
save(fullfile(savdir, 'hfb_ga_trl'), 'sleep', 'awake');

