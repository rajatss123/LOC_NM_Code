%% Format Acquired Data
load('C:\Users\NeuroLab\Documents\MATLAB\BITES\Subj_XX\acquired\CRAW_ecog.mat');
subj_char = 'Subj_XX'; 
acquire_dir = 'C:\Users\NeuroLab\Documents\MATLAB\BITES\';
dtype = 'ecog';
fs = 250;

pnt_start = 1;
pnts_ecog = locs_ecog*fs;

raw_data = {};
for idx = 1:length(pnts_ecog)
pnt_start = pnts_ecog(idx)+1;
pnt_end = pnt_start + 90*fs; %90s of data in each file
temp_data = cellfun(@(x) x(pnt_start:pnt_end), CRAW, 'UniformOutput', false);
raw_data{idx} = cat(1, temp_data{:});
end


ftdata = [];

ftdata.trial = raw_data;
ftdata.label = {'LeftNAc1-3' 'LeftNAc2-4' 'RightNAc1-3' 'RightNAc2-4'};
ftdata.fsample = 250;
ftdata.time = [0:1/ftdata.fsample:((size(ftdata.trial{1}, 2)-1)/ftdata.fsample)];
ftdata.time = {ftdata.time - 60};
ftdata.time = repelem(ftdata.time, length(ftdata.trial));

%% Process data

cfg = [];
cfg.lpfilter = 'yes';
cfg.lpfreq = 100;
cfg.demean = 'yes';
cfg.detrend = 'yes';

ftdata = ft_preprocessing(cfg, ftdata);

%%Epoch in 1 min pre-magnet swipe
cfg = [];
cfg.trials = 'all';
cfg.toilim = [-60 0];
ftdata = ft_redefinetrial(cfg,ftdata);

out_dir = [acquire_dir subj_char '\acquired_formated\pre60s_psd\'];

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

save([out_dir 'data_' dtype], '-struct', 'ftdata', '-v7.3');

%% Analyze PSD
% p = (gcp('nocreate'));
% 
% if isempty(p)
% parpool;
% end

%%Analyze data
cfg = [];
cfg.output = 'pow';
cfg.keeptrials = 'yes';
cfg.method = 'mtmfft';
cfg.foi = 1:0.25:90; 
cfg.taper = 'hanning'; 
cfg.pad = ceil(max(cellfun(@numel, ftdata.time)/ftdata.fsample));

ftdata_mfft_psd= ft_freqanalysis(cfg, ftdata);

out_dir = [subj_char '\analyzed_data\'];

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

save([out_dir 'pre60s_psd_' dtype], 'ftdata_mfft_psd', '-v7.3');

%% Calculate Mean and SE for each state
highcraving_trl = A:B; % Define trial(files) with High Craving, no hunger (LOC) trials. Replace letter with file idx number
hunger_trl = C:D; % Define trial(files) with Hunger. Replace letter with file idx number
normal_trl = E:F; % Replace letter with file idx number


mean_aggr = {};
sem_aggr = {};

mean_aggr{1} = squeeze(nanmean(ftdata_mfft_psd.powspctrm(hunger_trl, :, :), 1));
mean_aggr{2} = squeeze(nanmean(ftdata_mfft_psd.powspctrm(highcraving_trl, :, :), 1));
mean_aggr{3} = squeeze(nanmean(ftdata_mfft_psd.powspctrm(normal_trl, :, :),1));


sem_aggr{1} = squeeze(sem2(ftdata_mfft_psd.powspctrm(craving_hunger_trl, :, :), 1));
sem_aggr{2} = squeeze(sem2(ftdata_mfft_psd.powspctrm(highcraving_trl, :, :), 1));
sem_aggr{3} = squeeze(sem2(ftdata_mfft_psd.powspctrm(normal_trl, :, :),1));


save([out_dir 'pre60s_psd_mean_sem_' dtype], 'mean_aggr', 'sem_aggr', '-v7.3');

%% Plot PSD with SE error bar 
dir_acquired = 'C:\Users\NeuroLab\Documents\MATLAB\BITES\Subj_XX\';
date_now = datestr(datetime('now'), 'yyyy_mm_dd_HH_MM_SS');
out_dir = [dir_acquired 'plot_figures\' char(date_now) '\'];

    if ~exist(dir_acquired, 'dir')
        mkdir(dir_acquired);
    end
    
    if ~exist(out_dir)
        mkdir(out_dir);
    end


x_data = 1:0.25:90;
line_cols = {'b' 'r' 'k'}; 

%% Plot PSD w/SE for all 3 conditions

for i = 1:length(e_type)
    fig(i) = figure('units','normalized','outerposition',[0 0 1 1], 'Name', [subj_char '_' e_type{i}]);
    
    ck_shadedErrorBar(x_data, mean_aggr{1,1}(i,:), sem_aggr{1,1}(i,:), line_cols{1},1);
    hold on
    ck_shadedErrorBar(x_data, mean_aggr{1,2}(i,:), sem_aggr{1,2}(i,:), line_cols{2},1);
    ck_shadedErrorBar(x_data, mean_aggr{1,3}(i,:), sem_aggr{1,3}(i,:), line_cols{3},1); 
    
    xlabel('Frequency (Hz)');
    ylabel('Power (V^2)');
    xlim([1 90])
    title(e_type{i});
    hold off
    
    saveas(fig(i), [out_dir 'pre60s_PSD_allconditions_1to90Hz' char(strtrim(get(fig(i), 'name')))], 'png');
end
close all

for i = 1:length(e_type)
    fig(i) = figure('units','normalized','outerposition',[0 0 1 1], 'Name', [subj_char '_' e_type{i}]);
    
    ck_shadedErrorBar(x_data, mean_aggr{1,1}(i,:), sem_aggr{1,1}(i,:), line_cols{1},1);
    hold on
    ck_shadedErrorBar(x_data, mean_aggr{1,2}(i,:), sem_aggr{1,2}(i,:), line_cols{2},1);
    ck_shadedErrorBar(x_data, mean_aggr{1,3}(i,:), sem_aggr{1,3}(i,:), line_cols{3},1); 
    xlabel('Frequency (Hz)');
    ylabel('Power (V^2)');
    xlim([1 15])
    title(e_type{i});
    hold off
    
    saveas(fig(i), [out_dir 'pre60s_PSD_allconditions_1to15Hz' char(strtrim(get(fig(i), 'name')))], 'png');
end

        
                



        
        