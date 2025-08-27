function hfb_stat_lme_trl_reliability_sbj(num_sbj, num_trl)
% HFB_STAT_LME_TRL_RELIABILITY_SBJ - run linear mixed-effects model on HFB 
% sleep vs. wake data, varying the number subjects and the number of trials
% per subject, and repeating the analysis 10 times. Plot results on the 
% scalp topography. Outputs are FDR-corrected across channels using the 
% subfunction fdr.m.
%
% Ensure FieldTrip is correcty added to the MATLAB path:
%   addpath <path to fieldtrip home directory>
%   ft_defaults
%
% Inputs:
% num_sbj = number of subjects (e.g., 10)
% num_trl = number of trials (e.g., 5)
%
% Example:
% hfb_stat_lme_trl_reliability_sbj(10, 5)
%
% Copyright (c) 2023
% EL Johnson, PhD & AM Holubecki, PhD candidate

% set directories
pth = pwd;
datdir = fullfile(pth, 'All_Data_HFB'); % ouput of hfb_ga
savdir = fullfile(datdir, 'reliability', 'group_n');
mkdir(savdir);

num_x = 10; % # replications

% load data
load(fullfile(datdir, 'hfb_ga_trl'), 'sleep', 'awake');

% loop through replications
x_rand_s = cell(length(sleep.powspctrm), num_x);
x_rand_a = x_rand_s;
not_nan = zeros(length(sleep.powspctrm), length(sleep.label));
for x = 1:num_x
    for s = 1:length(sleep.powspctrm)
        % randomly sample trials without replacement
        x_rand_s{s,x} = randsample(size(sleep.powspctrm{s},1), num_trl, false);
        x_rand_a{s,x} = randsample(size(awake.powspctrm{s},1), num_trl, false);
        
        % index data channels per subject
        if x == 1
            for e = 1:length(sleep.label)
                if ~isnan(sleep.powspctrm{s}(1,e)) && ~isnan(awake.powspctrm{s}(1,e))
                    not_nan(s,e) = 1;
                end
            end
        end
    end
    not_nan = logical(not_nan);
    
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
        % confirm that the channel contains data
        if sum(not_nan(:,e)) >= num_sbj
            % extract data for table
            sbj = randsample(find(not_nan(:,e)), num_sbj, false);
            sid = [];
            trl = [];
            sleep_wake = [];
            dat = [];
            for s = 1:num_sbj
                tmp_dat = [sleep.powspctrm{sbj(s)}(x_rand_s{sbj(s),x},e); ...
                    awake.powspctrm{sbj(s)}(x_rand_a{sbj(s),x},e)];
                tmp_s_a = [zeros(size(sleep.powspctrm{sbj(s)}(x_rand_s{sbj(s),x},e))); ...
                    ones(size(awake.powspctrm{sbj(s)}(x_rand_a{sbj(s),x},e)))];
                tmp_sid = repmat(sbj(s), size(tmp_dat));
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
    end

    % apply FDR correction
    [~, lme.fdr_th] = fdr(lme.p, 0.05);
    lme.fdr = lme.p < lme.fdr_th;
    
    % save
    save(fullfile(savdir, ['hfb_stat_lme_n' num2str(num_sbj) '_' num2str(num_trl) 'trl_r' num2str(x)]), 'lme');

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
    print(fullfile(savdir, ['hfb_stat_lme_n' num2str(num_sbj) '_' num2str(num_trl) 'trl_r' num2str(x)]), '-dtiff');
end

end
