
load('C:\Users\NeuroLab\Documents\MATLAB\BITES_XX_BUFFET\Subj_XX_MMDDYYYY_Buffet\acquired\CRAW_ecog.mat');
subj_char = 'Subj_XX_MMDDYYYY_Buffet'; %XX update subject number, and date of task
dtype = 'ecog';  
acquire_dir = 'C:\Users\NeuroLab\Documents\MATLAB\BITES_XX_BUFFET\'; 

% % Add baseline data in front of trials 
% pnt_start = locs_ecog(length(locs_ecog))*fs;
% pnt_end = pnt_start + 10*fs;
% bsln_data = cellfun(@(x) x(pnt_start:pnt_end), CRAW, 'UniformOutput', false);

fs = 250;
pnt_start = 1;
pnts_ecog = locs_ecog*fs;

raw_data = {};
for idx = 1:length(pnts_ecog)
pnt_start = pnts_ecog(idx)+1;
pnt_end = pnt_start + 30*fs;
temp_data = cellfun(@(x) x(pnt_start:pnt_end), CRAW, 'UniformOutput', false);
raw_data{idx} = cat(1, temp_data{:});
end


ftdata = [];

ftdata.trial = raw_data;
ftdata.label = {'LNAc1-3' 'LNAc2-4' 'RNAc1-3' 'RNAc2-4'};
ftdata.fsample = 250;
ftdata.time = [0:1/ftdata.fsample:((size(ftdata.trial{1}, 2)-1)/ftdata.fsample)];
ftdata.time = {ftdata.time - 20};
ftdata.time = repelem(ftdata.time, length(ftdata.trial));

%% Pre-process data

cfg = [];
cfg.lpfilter = 'yes';
cfg.lpfreq = 100;
cfg.demean = 'yes';
cfg.detrend = 'yes';

ftdata = ft_preprocessing(cfg, ftdata);
ftdata.trialinfo = {};% Define Trials here LOC Bites, Normal Caloric Bites

out_dir = [acquire_dir subj_char '\acquired_formated\'];


if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

save([out_dir 'data_' dtype], '-struct', 'ftdata', '-v7.3');


%% Analyze PSD
cfg = [];
cfg.method = 'wavelet'; 
cfg.output = 'pow';
cfg.keeptrials = 'yes';
cfg.pad = 'nextpow2';
cfg.polyremoval = -1;


cfg.foi = 1:0.5:90;
cfg.taper = 'hanning'; 
cfg.t_ftimwin = ones(length(cfg.foi),1).*1; 
cfg.tapsmofrq = 6;
toi = [-20 10];
time_invl = .05;
time_sel = ftdata.time{1} >= toi(1) & ftdata.time{1} <= toi(2);
time_temp = ftdata.time{1}(time_sel);
cfg.toi = time_temp(1:ceil(time_invl*ftdata.fsample):length(time_temp));

ftdata_wavelet_psd= ft_freqanalysis(cfg, ftdata);

out_dir = [acquire_dir subj_char '\analyzed_data\'];

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

save([out_dir 'psd_' dtype], 'ftdata_wavelet_psd', '-v7.3');

%% Calculate 2D TF Spectral Plot
date_now = datestr(datetime('now'), 'yyyy_mm_dd_HH_MM_SS');
out_dir = [acquire_dir subj_char '\plot_figures\' char(date_now) '\'];

    if ~exist(out_dir)
        mkdir(out_dir);
    end

f_pref = 'TimeFrequencyData_lfp';
load([acquire_dir subj_char '\analyzed_data\test_psd_' dtype]);

for chan_idx = 1:length(ftdata_wavelet_psd.label)
    chan_char = ftdata_wavelet_psd.label{chan_idx};
    
    normal_trl = A:B; %Normal Meal Bites (replace letters with trial #)
    loc_trl = C:D; %LOC Meal Bites
    
    x_data = ftdata_wavelet_psd.time;
    y_data = ftdata_wavelet_psd.freq;
    
    fdata_aggr{1} = permute(nanmean(ftdata_wavelet_psd.powspctrm(normal_trl, chan_idx, :, :), 1), [3 4 1 2]);
    fdata_aggr{2} = permute(nanmean(ftdata_wavelet_psd.powspctrm(loc_trl, chan_idx, :, :), 1), [3 4 1 2]);
   

    % Relative Change
    fdata_aggr{3} = (fdata_aggr{1} - fdata_aggr{2})./fdata_aggr{1};
  
    
    f_type = {'Standard Meal Bites', 'HF Food - LOC',...
        'Standard vs HF Food Difference',};
    
    fig(chan_idx) = figure('units','normalized','outerposition',[0 0 1 1], 'Name', [subj_char chan_char]);
        plot_size = length(f_type)-2;
        plot_num = 1;
        for idx = 1:3
            subplot(plot_size,1, plot_num);
            plot_num = plot_num + 1;
            imagesc(x_data, y_data, fdata_aggr{idx}); 
            colormap(jet);
            hold on
            set(gca,'YDir','normal');
            colorbar;
            title(f_type{idx});
            hold off
        end
        saveas(fig(chan_idx), [out_dir 'TF_1-90Hz' char(strtrim(get(fig(chan_idx), 'name')))], 'png');
end


close all

%% Analyze TF
cfg = [];
cfg.method = 'mtmconvol'; 
cfg.output = 'pow';
cfg.keeptrials = 'yes';
cfg.foi = 1:1:90;
cfg.taper = 'dpss';
cfg.t_ftimwin = ones(length(cfg.foi),1).*0.5;
cfg.keeptapers = 'no';
cfg.tapsmofrq = 6; 
toi = [-20 10];
time_invl = .05;
time_sel = ftdata.time{1} >= toi(1) & ftdata.time{1} <= toi(2);
time_temp = ftdata.time{1}(time_sel);
cfg.toi = time_temp(1:ceil(time_invl*ftdata.fsample):length(time_temp));

ftdata_TF= ft_freqanalysis(cfg, ftdata);
out_dir = [acquire_dir subj_char '\analyzed_data\'];
save([out_dir 'test_tf_' dtype], 'ftdata_TF', '-v7.3');


%% PLOT TF Bandpower
  for chan_idx = 1:length(ftdata_bp.label)
        
        chan_char = ftdata_bp.label{chan_idx};
        
        cfg = [];
        cfg.channel = chan_char;
        outdata = ft_selectdata(cfg, ftdata_bp); %outdata

        
        cfg.baselinetype = 'zscore';
        cfg.baseline = [9 10];
        temp_data = ft_freqbaseline(cfg, outdata);
        
        cfg1 = [];
        plt_win = [-2 2]; %analyze 2s prior to 2s post bite. Window can vary based on video of bite onset
        cfg1.latency = plt_win;
        temp_data = ft_selectdata(cfg1, temp_data);
        
        normal_trl; 
        loc_trl; 
        
  
        
        x_data = temp_data.time;
        y_data = temp_data.freq;
        
        fdata_aggr ={};
        
        fdata_aggr{1} = permute(nanmean(temp_data.powspctrm(loc_trl, 1, :, :), 1), [3 4 1 2]);
        fdata_aggr{2} = permute(nanmean(temp_data.powspctrm(normal_trl, 1, :, :), 1), [3 4 1 2]);
        fdata_aggr{3} = fdata_aggr{1} - fdata_aggr{2};
        
        f_opt = [1 2; 2 8; 9 14; 15 39; 40 57; 57 90];
        
        f_label = {};
        mean_aggr = {};
        ste_aggr = {};
        for idx_temp = 1:size(f_opt, 1)
            f_label{idx_temp} = [char(num2str(f_opt(idx_temp, 1))) '-' char(num2str(f_opt(idx_temp, 2)))];
            f_sel = ftdata_bp.freq >= f_opt(idx_temp, 1) & ftdata_bp.freq <= f_opt(idx_temp, 2);
            
            data_aggr{1, idx_temp} = permute(mean(temp_data.powspctrm(loc_trl, 1, f_sel, :), 3), [1 4 2 3]);
            data_aggr{2, idx_temp} = permute(mean(temp_data.powspctrm(normal_trl, 1, f_sel, :), 3), [1 4 2 3]);
            
            
            
            
            for idx = 1:size(data_aggr, 1)
                mean_aggr{idx, idx_temp}  = mean(data_aggr{idx, idx_temp},1);
                ste_aggr{idx, idx_temp}  = sem2(data_aggr{idx, idx_temp},1);
            end
            
        end
        
        e_type = {'LOC' 'Normal' 'Difference'};
        plot_size = round(sqrt(3 + prod([size(mean_aggr, 2), 2])));
        set(0,'DefaultFigureVisible','off');
        fig = figure('units','normalized','outerposition',[0 0 1 1], 'Name', [subj_char chan_char]);
        plot_num = 1;
        for idx = 1:length(fdata_aggr)
            subplot(plot_size, plot_size+1, plot_num);
            plot_num = plot_num + 1;
            
            imagesc(x_data, y_data, fdata_aggr{idx}); colormap(jet); hold on
            set(gca,'YDir','normal');
            colorbar;
            ylim([1 90]);
            caxis([-10 10]);
            title(e_type{idx});
        end
        
        for freq_idx = 1:size(mean_aggr, 2)
            subplot(plot_size, plot_size+1, plot_num);
            plot_num = plot_num + 1;
            
            line_cols = {'r' 'k' 'g'};
            for idx = 1:size(mean_aggr, 1)
                ck_shadedErrorBar(x_data, mean_aggr{idx, freq_idx}, ste_aggr{idx, freq_idx}, line_cols{idx},1); hold on; box off;
            end
            
            y_max = ylim;
            poi = [0];
            
            for poi_idx = 1:length(poi)
                line([poi(poi_idx) poi(poi_idx)], y_max);
            end
            title(f_label{freq_idx});
            
            subplot(plot_size, plot_size+1, plot_num);
            plot_num = plot_num + 1;
            
            for idx = 1:size(data_aggr, 1)
                test_sigALL
            end
            
            test_sigDiffAll
        end
        

        set(fig,'PaperOrientation','portrait');
        saveas(fig, [out_dir 'Freq_' char(strtrim(get(fig, 'name')))], 'png');
        close all
        
  end
    






