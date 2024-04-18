function results_table = compare_psd(psd_struct1,psd_struct2)
% function to compare 2 psd structure
%
% Inputs :
%           - psd_struct : psd structure composed of N struct (N channels)
%               * channel
%               * fft_iterations
%               * psd
%               * f_array
%
%
% SMV : 06/03/2023

% %% inputs for debug : to comment!!!
% psd_struct1 = load('sub-663617f7d722_ses-c1_task-restEC_acq-j1_run-01_proc-fullcleaned_BF_on_VEPregions_fftwindow-256_fftoverlap-128_window-None.mat')
% f_array = 1: length(psd_struct1.psd(1).psd);
% psd_struct2 = load('sub-663617f7d722_ses-c1_task-restEC_acq-j1_run-02_proc-fullcleaned_BF_on_VEPregions_fftwindow-256_fftoverlap-128_window-None.mat')
% psd_struct1 = psd_struct1.psd;
% for i = 1 : length(psd_struct1)
%     psd_struct1(i).f_array = f_array;
% end
% psd_struct2 = psd_struct2.psd;
% for i = 1 : length(psd_struct2)
%     psd_struct2(i).f_array = f_array;
% end

%% check inputs :
channels1 = {psd_struct1(:).channel};
channels2 = {psd_struct2(:).channel};

chan_to_do = intersect(channels1,channels2);

% f_band = [1 4; 4 8; 8 15; 15 30; 30 45; 45 80; 80 250];
f_band = [4 8; 8 15; 15 30; 30 45];

[M,f4_idx1 ]= min(abs(psd_struct1(1).f_array-f_band(1,1)));
[M,f45_idx2 ]= min(abs(psd_struct1(1).f_array-f_band(4,2)));

results_table = struct;
a = 0;
%% Computation
nb_chan = length(chan_to_do);
mean_spectrum1 = [];
mean_spectrum2 = [];
percentage_band_power_tot1 =[];
percentage_band_power_tot2 = [];
fooof_band_power_tot1 = [];
fooof_band_power_tot2 = [];
for i = 1 : nb_chan
    idx1 = find(strcmpi({psd_struct1.channel},chan_to_do{i}))
    idx2 = find(strcmpi({psd_struct2.channel},chan_to_do{i}))
    %     %% display spectrum of chan 1 et 2
%         f1 = figure;
%         hold on;
%         plot(psd_struct1(idx1).f_array,psd_struct1(idx1).psd)
%         plot(psd_struct2(idx2).f_array,psd_struct2(idx2).psd)
%         legend({'stim','baseline'})
%         title(chan_to_do{i},'Interpreter','none')
%         print(f1, ['figure ' chan_to_do{i}], '-djpeg')
%         close(f1)

    %
%         figure;
%         hold on;
%         plot(psd_struct1(idx1).f_array,psd_struct1(idx1).psd)
%        mean_psd = mean(psd_struct1(idx1).fft_iterations(:,:),2)
%         plot(psd_struct1(idx1).f_array,10*log10(mean_psd))
    %% Compute on mean spectrum 
    % FOOOF
    % FOOOF settings
    freqs = psd_struct1(idx1).f_array;
    f_range = [f_band(1,1) f_band(4,2)];
    settings = struct();  % Use defaults
    freq_oi = freqs(f4_idx1:f45_idx2);
    idx_freq_oi = f4_idx1:f45_idx2;

    mean_spectrum1(:,i) =  mean(psd_struct1(idx1).fft_iterations,2);
    mean_spectrum2(:,i) =  mean(psd_struct2(idx2).fft_iterations,2);
    total_meanpow1 = sum(mean_spectrum1(f4_idx1:f45_idx2,i),1);
    total_meanpow2= sum(mean_spectrum2(f4_idx1:f45_idx2,i),1);
    
    if i == 1
        figure;plot([mean_spectrum1(:,1)/total_meanpow1,mean_spectrum2(:,1)/total_meanpow2])
        title([chan_to_do{i} ' mean spectrum'],'Interpreter','none')
    end

    oscill_part_tot1 = zeros(size(psd_struct1(idx1).f_array));
    oscill_part_tot2 = zeros(size(psd_struct1(idx1).f_array));
    % fooof ap
    fooof_results_tot1 = fooof(freqs, mean_spectrum1(:,i), f_range, settings,true);
    fooof_results_tot2 = fooof(freqs, mean_spectrum2(:,i), f_range, settings,true);

    if i ==1
        fooof_plot(fooof_results_tot1)
        title([chan_to_do{i} ' fooof plot'],'Interpreter','none')
    end


    oscill_part_tot1(idx_freq_oi) = fooof_results_tot1.fooofed_spectrum-fooof_results_tot1.ap_fit;
    aperiodic_param_tot1(i) = fooof_results_tot1.aperiodic_params(2);

    oscill_part_tot2(idx_freq_oi) = fooof_results_tot2.fooofed_spectrum-fooof_results_tot2.ap_fit;
    aperiodic_param_tot2(i) = fooof_results_tot2.aperiodic_params(2);

    if i ==1
        figure; plot(oscill_part_tot1); hold on; plot(oscill_part_tot2)
        title([chan_to_do{i} ' fooof oscillatory parts'],'Interpreter','none')
    end
    % SLOPE
    [b1,stats1] = robustfit(freqs(f4_idx1:f45_idx2),mean_spectrum1(f4_idx1:f45_idx2,i));
    [b2,stats2] = robustfit(freqs(f4_idx1:f45_idx2),mean_spectrum2(f4_idx1:f45_idx2,i));
    slope_tot1(i) = b1(2);
    slope_tot2(i) = b2(2);
    if i ==1
        figure; hold on;
        plot(freqs(f4_idx1:f45_idx2),mean_spectrum1(f4_idx1:f45_idx2,i))
        plot(freqs(f4_idx1:f45_idx2), b1(1) + b1(2)*freqs(f4_idx1:f45_idx2))
        plot(freqs(f4_idx1:f45_idx2),mean_spectrum2(f4_idx1:f45_idx2,i))
        plot(freqs(f4_idx1:f45_idx2), b2(1) + b2(2)*freqs(f4_idx1:f45_idx2))
        title([chan_to_do{i} ' mean spectrums and slopes'],'Interpreter', 'none')
    end
    %% Compare psd_iterations per band
    total_pow1 = sum(psd_struct1(idx1).fft_iterations(f4_idx1:f45_idx2,:),1);           % regarder entre 4 et 45
    total_pow2 = sum(psd_struct2(idx2).fft_iterations(f4_idx1:f45_idx2,:),1);
    percentage_band_power1 = [];
    percentage_band_power2 = [];
    %% faire fooof sur chaque fenetre
    aperiodic_param1 = [];
    aperiodic_param2 = [];
    nb_win1 = size(psd_struct1(idx1).fft_iterations,2);
    nb_win2 = size(psd_struct2(idx2).fft_iterations,2);
    oscill_part1 = zeros(size(psd_struct1(idx1).f_array,2),nb_win1);
    oscill_part2 = zeros(size(psd_struct1(idx1).f_array,2),nb_win2);
    parfor w = 1:nb_win1
        % Transpose, to make inputs row vectors
        psd = psd_struct1(idx1).fft_iterations(:,w)';

        %         % Run FOOOF
        %         fooof_results = fooof(freqs, psd_struct1(idx1).fft_iterations(:,1), f_range, settings,true);
        %         % Print out the FOOOF Results
        %         fooof_plot(fooof_results)
        % Run FOOOF
        fooof_results1 = fooof(freqs, psd, f_range, settings,true);

        oscill_part1(idx_freq_oi,w) = fooof_results1.fooofed_spectrum-fooof_results1.ap_fit;
        aperiodic_param1(w) = fooof_results1.aperiodic_params(2);
    end
    parfor w = 1:nb_win2
        % Transpose, to make inputs row vectors
        psd = psd_struct2(idx2).fft_iterations(:,w)';

        % Run FOOOF
        fooof_results2 = fooof(freqs, psd, f_range, settings,true);

        oscill_part2(idx_freq_oi,w) = fooof_results2.fooofed_spectrum-fooof_results2.ap_fit;
        aperiodic_param2(w) = fooof_results2.aperiodic_params(2);
    end
    % Make the comparison of energy per band
    percentage_band_power1 = [];
    percentage_band_power2 = [];
    fooof_band_power1 = [];
    fooof_band_power2 = [];

    for f_band_idx = 1 : size(f_band,1)
        a = a+1;
        [M,f_idx1 ]= min(abs(psd_struct1(idx1).f_array-f_band(f_band_idx,1)));
        [M,f_idx2 ]= min(abs(psd_struct1(idx1).f_array-f_band(f_band_idx,2)));

        %% compute fooof and band power on average spectrum
        % tot_band_power
        [M,f_idx1 ]= min(abs(freqs-f_band(f_band_idx,1)));
        [M,f_idx2 ]= min(abs(freqs-f_band(f_band_idx,2)));
        
        percentage_band_power_tot1(i,f_band_idx) = sum(mean_spectrum1(f_idx1:f_idx2,i),1)./total_meanpow1;
        percentage_band_power_tot2(i,f_band_idx) = sum(mean_spectrum2(f_idx1:f_idx2,i),1)./total_meanpow2;

        % fooof band power
        fooof_band_power_tot1(f_band_idx) = mean(oscill_part_tot1(f_idx1:f_idx2));
        fooof_band_power_tot2(f_band_idx) = mean(oscill_part_tot2(f_idx1:f_idx2));

        %% comparison of band power
        percentage_band_power1(f_band_idx,:) = sum(psd_struct1(idx1).fft_iterations(f_idx1:f_idx2,:),1)./total_pow1;
        percentage_band_power2(f_band_idx,:) = sum(psd_struct2(idx2).fft_iterations(f_idx1:f_idx2,:),1)./total_pow2;
        %         band_power1(f_band_idx,:) = mean(psd_struct1(idx1).fft_iterations(f_idx1:f_idx2,:),1);
        %         band_power2(f_band_idx,:) = mean(psd_struct2(idx2).fft_iterations(f_idx1:f_idx2,:),1);

        [p_band_pow_comp,h_band_pow_comp,stats_band_pow_comp] = ranksum(percentage_band_power1(f_band_idx,:)',percentage_band_power2(f_band_idx,:)');

        log_spectrum_by_win1 = 10*log10(psd_struct1(idx1).fft_iterations);
        log_spectrum_by_win2 = 10*log10(psd_struct2(idx2).fft_iterations);

        %visu
        if i == 1
            figure; histogram(percentage_band_power1(f_band_idx,:),50); hold on; ...
                histogram(percentage_band_power2(f_band_idx,:),50);
            title([chan_to_do{i} ['BANDPOW hist of band power band ' num2str(freqs(f_idx1)) ...
                ' but p='] num2str(p_band_pow_comp)],'Interpreter','none')
        end

        %% FOOOF oscillatory comparison
        fooof_band_power1(f_band_idx,:) = mean(oscill_part1(f_idx1:f_idx2,:),1);
        fooof_band_power2(f_band_idx,:) = mean(oscill_part2(f_idx1:f_idx2,:),1);
        [p_fooof_osc_comp,h_fooof_osc_comp,stats_fooof_osc_comp] = ranksum(fooof_band_power1(f_band_idx,:)',fooof_band_power2(f_band_idx,:)');
        
        %visu
        if i == 1
            figure; histogram(fooof_band_power1(f_band_idx,:),50); hold on; histogram(fooof_band_power2(f_band_idx,:),50);
            title([chan_to_do{i} ['FOOOF hist of band power band ' num2str(freqs(f_idx1)) ...
                ' but p='] num2str(p_fooof_osc_comp)],'Interpreter','none')
        end

        %% OUTPUT RESULTS
        results_table(a).channel = chan_to_do{i};
        results_table(a).band = f_band(f_band_idx,:);
        results_table(a).p_band_pow_comp = p_band_pow_comp;
        results_table(a).z_band_pow_comp= stats_band_pow_comp.zval;
        results_table(a).percentage_band_power1 = percentage_band_power1(f_band_idx,:);
        results_table(a).percentage_band_power2 = percentage_band_power2(f_band_idx,:);
        %         results_table(a).band_power1 = band_power1(f_band_idx,:);
        %         results_table(a).band_power2 = band_power2(f_band_idx,:);
        results_table(a).p_fooof_band_power = p_fooof_osc_comp;
        results_table(a).z_fooof_band_power = stats_fooof_osc_comp.zval;
        results_table(a).fooof_band_power1 = fooof_band_power1(f_band_idx,:);
        results_table(a).fooof_band_power2 = fooof_band_power2(f_band_idx,:);
        % write results for average analysis
        results_table(a).fooof_band_power_tot1 = fooof_band_power_tot1(f_band_idx);
        results_table(a).fooof_band_power_tot2 = fooof_band_power_tot2(f_band_idx);
        results_table(a).percentage_band_power_tot1 = percentage_band_power_tot1(i,f_band_idx);
        results_table(a).percentage_band_power_tot2 = percentage_band_power_tot2(i,f_band_idx);
    end
%% average spectrum results output
results_table(a).aperiodic_param_tot1 = aperiodic_param_tot1(i);
results_table(a).aperiodic_param_tot2 = aperiodic_param_tot2(i);
results_table(a).slope_tot1 = slope_tot1(i);
results_table(a).slope_tot2 = slope_tot2(i);

    %% fooof outputs
    % aperiodic comp 
%     figure;
%     hold on;
%     histogram(aperiodic_param1,30)
%     histogram(aperiodic_param2,30)
%     hold off;
%     title(chan_to_do{i},'Interpreter','none')
% figure; boxplot(aperiodic_param1);ylim([-2 5])
% figure; boxplot(aperiodic_param2);ylim([-2 5])

    [p,h,stats] = ranksum(aperiodic_param1,aperiodic_param2);
%     [h,p] = ttest2(aperiodic_param1,aperiodic_param2)

    results_table(a).aperiodic_comp_p = p;
    results_table(a).aperiodic_comp_z = stats.zval;

    %% faire la pente sam
    stats1 = {};
    stats2 = {};
    b1 = [];
    b2 = [];
    for w = 1:size(log_spectrum_by_win1,2)
        [b,stats] = robustfit(psd_struct1(idx1).f_array(f4_idx1:f45_idx2)',log_spectrum_by_win1(f4_idx1:f45_idx2,w));
        b1(w,:) = b;
        stats1(w) = {stats};
        all_slope1 = b1(:,2);
    end
    for w = 1:size(log_spectrum_by_win2,2)
        [b,stats] = robustfit(psd_struct2(idx2).f_array(f4_idx1:f45_idx2)',log_spectrum_by_win2(f4_idx1:f45_idx2,w));
        b2(w,:) = b;
        stats2(w) = {stats};
        all_slope2 = b2(:,2);
    end
    % faire la stat
    [p,h,stats] = ranksum(all_slope1,all_slope2);
    results_table(a).slopes_comp_p = p;
    results_table(a).slopes_comp_z = stats.zval;

    if i ==1
        figure;
        hold on;
        plot(log_spectrum_by_win1(f4_idx1:f45_idx2,50));
        plot(b1(50,1)+b1(50,2)*psd_struct1(idx1).f_array(f4_idx1:f45_idx2))
        plot(log_spectrum_by_win2(f4_idx1:f45_idx2,50));
        plot(b2(50,1)+b2(50,2)*psd_struct1(idx1).f_array(f4_idx1:f45_idx2))

        figure; 
        hold on;
        plot(mean(log_spectrum_by_win1(f4_idx1:f45_idx2,:),2));
        plot(mean(b1(:,1),1)+mean(b1(:,2),1)*psd_struct1(idx1).f_array(f4_idx1:f45_idx2))
        plot(mean(log_spectrum_by_win2(f4_idx1:f45_idx2,:),2));
        plot(mean(b2(:,1),1)+mean(b2(:,2),1)*psd_struct1(idx1).f_array(f4_idx1:f45_idx2))
        title([chan_to_do{i} ' average psd but p=' num2str(p)],'Interpreter','none')

        % visu
        figure;
        hold on;
        histogram(all_slope1,100)
        histogram(all_slope2,100)
        hold off;
        title([chan_to_do{i} ' average psd but p=' num2str(p)],'Interpreter','none')
    end

end



%% test sam slope accross all channels
% % find slope
% Y1  = psd_struct1(idx1).psd;
% [b1,stats1] = robustfit(1:length(Y1),Y1)
% Y2  = psd_struct2(idx2).psd;
% [b2,stats2] = robustfit(1:length(Y2),Y2)
% figure;
% hold on;
% plot(Y1);
% plot(b1(1)+b1(2)*[1:length(Y1)])
% plot(Y2);
% plot(b2(1)+b2(2)*[1:length(Y2)])
% 
% figure;
% s1 = subplot(2,1,1);
% s2 = subplot(2,1,2);
% for i = 1 : length(psd_struct1)
%     Y1  = psd_struct1(i).psd;
%     [b1,stats1] = robustfit(1:length(Y1),Y1);
%     all_slope(i,1) = b1(2);
%     Y2  = psd_struct2(i).psd;
%     [b2,stats2] = robustfit(1:length(Y2),Y2);
%     all_slope(i,2) = b2(2);
%     hold(s1, 'on');
%     plot(s1,b1(1)+b1(2)*[1:length(Y1)])
%     hold(s1, 'off');
%     hold(s2, 'on');
%     plot(s2,b2(1)+b2(2)*[1:length(Y2)])
%     hold(s2, 'off');
% end
% 
% % stat :
% [p,h,stats] = ranksum(all_slope(i,1),all_slope(i,2));
% % display
% figure;
% hold on;
% histogram(all_slope(:,1),30)
% histogram(all_slope(:,2),30)
% title(['histogram of all channels slope for 1 and 2. Wilcoxon p='...
%     num2str(p)])

end


