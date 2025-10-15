% HFB_ANALYSIS - compute high-frequency broadband (HFB; 70-150 Hz) power.
%
% Ensure FieldTrip is correcty added to the MATLAB path:
%   addpath <path to fieldtrip home directory>
%   ft_defaults
%
% Copyright (c) 2023
% EL Johnson, PhD & AM Holubecki, PhD candidate

% set directories
pth = pwd;
datdir = fullfile(pth, 'All_Data_clean_z');
savdir = fullfile(pth, 'All_Data_HFB');
mkdir(savdir);

% set up HFB filter
cfg = [];
cfg.pad = 'nextpow2';
cfg.method = 'mtmfft'; 
cfg.taper = 'dpss'; 
cfg.foi = 75:10:145;
cfg.tapsmofrq = 10/2;
cfg.output = 'pow'; % power
cfg.keeptrials = 'yes';

% loop through datasets
files = dir(fullfile(datdir, 'Subj*'));
for s = 1:length(files)
    % load data
    load(fullfile(datdir, files(s).name));

    % compute HFB
    sleep_eeg = ft_freqanalysis(cfg, sleep_eeg);
    awake_eeg = ft_freqanalysis(cfg, awake_eeg);

    % normalize each frequency band
    sleep_norm2 = nan(size(sleep_eeg.powspctrm));
    awake_norm2 = nan(size(awake_eeg.powspctrm));
    for f = 1:length(sleep_eeg.freq)
        sleep_norm2(:,:,f) = sleep_eeg.powspctrm(:,:,f) .* sleep_eeg.freq(f);
        awake_norm2(:,:,f) = awake_eeg.powspctrm(:,:,f) .* awake_eeg.freq(f);
    end

    % average across frequency bands
    sleep_eeg.powspctrm = mean(sleep_norm2,3);
    awake_eeg.powspctrm = mean(awake_norm2,3);

    sleep_eeg.freq = mean(sleep_eeg.freq);
    awake_eeg.freq = mean(awake_eeg.freq);

    % save
    save(fullfile(savdir, files(s).name), '*eeg');

    clear *norm2 *eeg
end


