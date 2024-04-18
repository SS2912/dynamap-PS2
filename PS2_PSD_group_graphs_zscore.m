%%%%%%  PSD PS2  %%%%%%%
% script to rin graphs for single subjects and group figures for all tDCS
% and all tACS subjects dividng into targeted and onn-targeted (same
% definition as in R) 

%% Graphs with zvalues (Pooled subjects)
opengl('save', 'software')
clear all
close all
format long g

% define parameters
n= 10;   % number of channels to average for each subj and each target
thrEF = 0.001;
tacs = struct();
tdcs = struct();
label = ["sham", "stimA", "stimB", "postA", "postB", "postALL"];
window   = 5; 
baseline = "baseline";
subbands = [1 4; 4 8; 8 15; 15 30; 30 45; 1 45]; 
subbands_name = ["delta", "theta", "alpha", "beta", "lowgamma", "broad"];
fmax   = max(max(subbands));


% upload table
data = load("C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\PSD\2023-11-15_13.53_PSDtableVSbaseline_fmax45.mat");
subject = data.PSDall;
tdcs_count = 0;
tacs_count = 0;
tdcs_stats_inhib = [];
tdcs_stats_exc = [];
tdcs_stats_na = [];
tacs_stats_targ = [];
tacs_stats_na = [];

% for each subject, select top n channels targeted or non targeted
for subj = 1:length(subject)
    tmp = subject(subj).chaninfo;
    tmp = tmp(strcmp(tmp.matter, "grey"),:);  % only take grey matter channels
    simEFnum = str2double(tmp.simEF);
    
    switch subject(subj).group

        case "tDCS"
        tdcs_count = tdcs_count +1;

        tdcs_exc = tmp(:, ["channel", "simEF"]);
        [~, idx] = sort(str2double(tdcs_exc.simEF), 'descend');
        tdcs_exc_chan = table2array(tdcs_exc(idx(1:n), "channel")); clearvars idx % take only the first n most inhib channels (to avoid bias in number of channels across subjs)
        idx_exc = ismember(subject(subj).chan, tdcs_exc_chan);

        tdcs_inhib = tmp(:, ["channel", "simEF"]);
        [~, idx] = sort(str2double(tdcs_inhib.simEF), 'ascend');
        tdcs_inhib_chan = table2array(tdcs_inhib(idx(1:n), "channel")); clearvars idx
        idx_inhib = ismember(subject(subj).chan, tdcs_inhib_chan);

        tdcs_na = tmp(:, ["channel", "simEF"]);
        [~, idx] = sort(abs(str2double(tdcs_na.simEF)), 'ascend');
        tdcs_na_chan = table2array(tdcs_na(idx(1:n), "channel")); clearvars idx
        idx_tdcsna = ismember(subject(subj).chan, tdcs_na_chan);

        % select data to plot from big structure
        tdcs(tdcs_count).subj = subject(subj).code;
        tdcs(tdcs_count).group = subject(subj).group;
        tdcs(tdcs_count).psd_inhib = {subject(subj).zvalSHAM(:,idx_inhib)', subject(subj).zvalSTIMA(:,idx_inhib)', subject(subj).zvalSTIMB(:,idx_inhib)', ...
            subject(subj).zvalPOSTA(:,idx_inhib)', subject(subj).zvalPOSTB(:,idx_inhib)', subject(subj).zvalPOSTALL(:,idx_inhib)'};
        tdcs(tdcs_count).psd_exc = {subject(subj).zvalSHAM(:,idx_exc)', subject(subj).zvalSTIMA(:,idx_exc)', subject(subj).zvalSTIMB(:,idx_exc)', ...
            subject(subj).zvalPOSTA(:,idx_exc)', subject(subj).zvalPOSTB(:,idx_exc)', subject(subj).zvalPOSTALL(:,idx_exc)'};
        tdcs(tdcs_count).psd_na = {subject(subj).zvalSHAM(:,idx_tdcsna)', subject(subj).zvalSTIMA(:,idx_tdcsna)', subject(subj).zvalSTIMB(:,idx_tdcsna)', ...
            subject(subj).zvalPOSTA(:,idx_tdcsna)', subject(subj).zvalPOSTB(:,idx_tdcsna)', subject(subj).zvalPOSTALL(:,idx_tdcsna)'};
        
        % build a 4D matrix with rows= mean of psd by subject across selected channels (n), cols=freqs, z=inhib or exc
        % or na and 4th-d=section (6 total)
        for section = 1:length(label)
%             mean_tdcs(tdcs_count, :, 1, section) = mean(tdcs(tdcs_count).psd_inhib{1,section}); 
%             mean_tdcs(tdcs_count, :, 2, section) = mean(tdcs(tdcs_count).psd_exc{1,section});
%             mean_tdcs(tdcs_count, :, 3, section) = mean(tdcs(tdcs_count).psd_na{1,section});
            mean_tdcs(tdcs_count, :, 1, section) = median(tdcs(tdcs_count).psd_inhib{1,section}); 
            mean_tdcs(tdcs_count, :, 2, section) = median(tdcs(tdcs_count).psd_exc{1,section});
            mean_tdcs(tdcs_count, :, 3, section) = median(tdcs(tdcs_count).psd_na{1,section});

            tmp_inhib(:,:,section)  = tdcs(tdcs_count).psd_inhib{1,section};
            tmp_exc(:,:,section)  = tdcs(tdcs_count).psd_exc{1,section};
            tmp_na(:,:,section)  = tdcs(tdcs_count).psd_na{1,section};
        end

        tdcs_stats_inhib(tdcs_count, :, :, :) = tmp_inhib;
        tdcs_stats_exc(tdcs_count, :, :, :) = tmp_exc;
        tdcs_stats_na(tdcs_count, :, :, :) = tmp_na;
       
        case "tACS"
        tacs_count = tacs_count+1;

        tacs_targ = tmp(:, ["channel", "simEF"]);
        [~, idx] = sort(abs(str2double(tacs_targ.simEF)), 'descend');
        tacs_targ_chan = table2array(tacs_targ(idx(1:n), "channel")); clearvars idx
        idx_targ = ismember(subject(subj).chan, tacs_targ_chan);
    
        tacs_na = tmp(:, ["channel", "simEF"]);
        [~, idx] = sort(abs(str2double(tacs_na.simEF)), 'ascend');
        tacs_na_chan = table2array(tacs_na(idx(1:n), "channel")); clearvars idx
        idx_tana = ismember(subject(subj).chan, tacs_na_chan);

        tacs(tacs_count).subj = subject(subj).code;
        tacs(tacs_count).group = subject(subj).group;
        tacs(tacs_count).psd_targ  = {subject(subj).zvalSHAM(:,idx_targ)', subject(subj).zvalSTIMA(:,idx_targ)', subject(subj).zvalSTIMB(:,idx_targ)', ...
            subject(subj).zvalPOSTA(:,idx_targ)', subject(subj).zvalPOSTB(:,idx_targ)', subject(subj).zvalPOSTALL(:,idx_targ)'};
        tacs(tacs_count).psd_na = {subject(subj).zvalSHAM(:,idx_tana)', subject(subj).zvalSTIMA(:,idx_tana)', subject(subj).zvalSTIMB(:,idx_tana)', ...
            subject(subj).zvalPOSTA(:,idx_tana)', subject(subj).zvalPOSTB(:,idx_tana)', subject(subj).zvalPOSTALL(:,idx_tana)'};

        for section = 1:length(label)
%             mean_tacs(tacs_count, :, 1, section) = mean(tacs(tacs_count).psd_targ{1,section}); 
%             mean_tacs(tacs_count, :, 2, section) = mean(tacs(tacs_count).psd_na{1,section});
            mean_tacs(tacs_count, :, 1, section) = median(tacs(tacs_count).psd_targ{1,section}); 
            mean_tacs(tacs_count, :, 2, section) = median(tacs(tacs_count).psd_na{1,section});

            tmp_targ(:,:,section) = tacs(tacs_count).psd_targ{1,section};
            tmp_tana(:,:,section) = tacs(tacs_count).psd_na{1,section};
        end

        tacs_stats_targ(tacs_count, :, :, :) = tmp_targ;
        tacs_stats_na(tacs_count, :, :, :) = tmp_tana;
                   

    end

end

%% stats z-values compared to zero
bands    = [1 4; 4 8; 8 13; 13 30; 30 45]; % ONLY WORKS WITH 45 max now, otherwise change it in hard code
band_name = ["delta", "theta", "alpha", "beta", "gamma"];

data_all = {tdcs_stats_inhib,tdcs_stats_exc, tdcs_stats_na, tacs_stats_targ, tacs_stats_na};
for i = 1:length(data_all)
    cur_data = data_all{i};
    for ii = 1:size(data_all{i},4)  % 4th dim = sections
%         median_cur_data = median(cur_data,2);  % median across channels for each subj
%         median_cur_data = squeeze(median_cur_data);
        for j=1:length(bands)
            f_min = bands(j,1)*window;
            f_max = bands(j,2)*window;
            cur_band = squeeze(median(cur_data(:,:,f_min:f_max,ii), 2)); 
            for sub = 1:size(cur_band,1)
                p(sub,j) = signtest(cur_band(sub,:)');
            end
        end
%             p(i,j,ii) = anova1(median_cur_data(:,:,ii)'); %anova on subjects (columns), data=median zvalue PSD, all frequencies pooled (rows)
        stats(ii).section = label(ii);
        stats(ii).data = i;
        stats(ii).pval = p;
        clearvars p % TO CHECK IF THIS WORKSS
    end
end


%%
supermean_tdcs = squeeze(median(mean_tdcs));
supermean_tacs = squeeze(median(mean_tacs)); 
colorlim = [-4 4];
color = viridis;


close all

figure('Name', strcat('tDCS: z-values compared with ', baseline, '-', int2str(n), 'channels median'))
for section = 1:length(label)
    % stats on the mean for pooled subjects
    disp(label(section))
    [p_inhib, h_inh, stats_inhib] = signtest(supermean_tdcs(:,1,section));
    [p_exc, h_exc, stats_exc] = signtest(supermean_tdcs(:,2,section));
    [p_na, h_na, stats_na] = signtest(supermean_tdcs(:,3,section));
    subplot(6,1,section)
    imagesc(supermean_tdcs(:,:,section)')
    colormap(color)
    xticks(0:fmax*window/10:fmax*window);
    xticklabels(0:fmax/10:fmax);
    yticks(1:n*3);
    yticklabels(["inhib", "exc", "na"]);
    colorbar
    clim(colorlim)
    if strcmp(label(section), "sham") && strcmp(baseline, "sham")
        title(strcat(label(section), " vs baseline"));
    else 
        title(strcat(label(section), " vs ", baseline));
    end
    xlab = "Frequency (Hz)";
    ylab = "Channel target type";
    ylabel(ylab);
    xlabel(xlab);
end

%
% tACS
figure('Name', strcat('tACS: z-values compared with ', baseline, '-', int2str(n), 'channels median'))

for section = 1:length(label)
    subplot(6,1,section)
    imagesc(supermean_tacs(:,:,section)')
    colormap(color)
    xticks(0:fmax*window/10:fmax*window);
    xticklabels(0:fmax/10:fmax);
    yticks(1:n);
    yticklabels(["targeted", "na"]);
    colorbar
    clim(colorlim)
    if strcmp(label(section), "sham") && strcmp(baseline, "sham")
        title(strcat(label(section), " vs baseline"));
    else 
        title(strcat(label(section), " vs ", baseline));
    end
    xlab = "Frequency (Hz)";
    ylab = "Channel target type";
    ylabel(ylab);
    xlabel(xlab);
end



%% 

close all   
clearvars subj

for subj = 1:length(tdcs)
figure('Name', strcat(tdcs(subj).subj, '(', tdcs(subj).group, ')', ...
        ': z-values compared with ', baseline))
    for section = 1:length(label)

        data = [tdcs(subj).psd_inhib{1,section}; tdcs(subj).psd_exc{1,section}; tdcs(subj).psd_na{1,section}]';
        
        subplot(2,3,section)
        imagesc(data)
        colormap(color)
        yticks(0:fmax*window/10:fmax*window);
        yticklabels(0:fmax/10:fmax);
        xticks(1:n*3);
        xticklabels(repelem(["inhib", "exc", "na"], n));
        colorbar
        clim(colorlim)
        if strcmp(label(section), "sham") && strcmp(baseline, "sham")
            title(strcat(label(section), " vs baseline"));
        else 
            title(strcat(label(section), " vs ", baseline));
        end
        ylabel(ylab);
        xlabel(xlab);
    end
end
    

    


%% Graphs with zvalues (1 figure per patient)
opengl('save', 'software')

close all
colorlim = [-10 10];
% when you need to charge older analysis:
% subject = load("C:\Users\simula\OneDrive - Aix-Marseille
% Université\PhD\MyProjects\Galvani\analysis\PSD\2023-04-27_18.28_PSDtableVSsham.mat");


for subj = 1:length(subject)

    data  = {subject(subj).zvalSHAM', subject(subj).zvalSTIMA', subject(subj).zvalSTIMB', ...
        subject(subj).zvalPOSTA', subject(subj).zvalPOSTB', subject(subj).zvalPOSTALL'};
    label = ["sham", "stimA", "stimB", "postA", "postB", "postALL"];
    figure('Name', strcat(subject(subj).code, '(', subject(subj).group, ')', ...
        ': z-values compared with ', baseline))
    
    for i=1:length(data)
        subplot(2,3,i)
        imagesc(data{i})
        colormap(color)
        xticks(0:fmax*window/10:fmax*window);
        xticklabels(0:fmax/10:fmax);
        yticks(1:length(subject(subj).chan'));
        yticklabels(subject(subj).chan');
        colorbar
        clim(colorlim)
        if strcmp(label(i), "sham") && strcmp(baseline, "sham")
            title(strcat(label(i), " vs baseline"));
        else 
            title(strcat(label(i), " vs ", baseline));
        end
        ylabel('channel label');
        xlabel('frequency (Hz)');
    end

end




%%%%% OLD
% tdcs_exc = tmp(str2double(tmp.simEF) > thrEF, ["channel", "simEF"]);
%         [~, idx] = sort(str2double(tdcs_exc.simEF), 'descend');
%         tdcs_exc_chan = table2array(tdcs_exc(idx(1:n), "channel")); clearvars idx % take only the first n most inhib channels (to avoid bias in number of channels across subjs)
%         idx_exc = ismember(subject(subj).chan, tdcs_exc_chan);
% 
%         tdcs_inhib = tmp(str2double(tmp.simEF) < 0 - thrEF, ["channel", "simEF"]);
%         [~, idx] = sort(str2double(tdcs_inhib.simEF), 'descend');
%         tdcs_inhib_chan = table2array(tdcs_inhib(idx(1:n), "channel")); clearvars idx
%         idx_inhib = ismember(subject(subj).chan, tdcs_inhib_chan);
% 
%         tdcs_na = tmp(abs(str2double(tmp.simEF)) < thrEF, ["channel", "simEF"]);
%         [~, idx] = sort(abs(str2double(tdcs_na.simEF)), 'ascend');
%         tdcs_na_chan = table2array(tdcs_na(idx(1:n), "channel")); clearvars idx
%         idx_tdcsna = ismember(subject(subj).chan, tdcs_na_chan);
% 
%         % select data to plot from big structure
%         tdcs(tdcs_count).subj = subject(subj).code;
%         tdcs(tdcs_count).group = subject(subj).group;
%         tdcs(tdcs_count).psd_inhib = {subject(subj).zvalSHAM(:,idx_inhib)', subject(subj).zvalSTIMA(:,idx_inhib)', subject(subj).zvalSTIMB(:,idx_inhib)', ...
%             subject(subj).zvalPOSTA(:,idx_inhib)', subject(subj).zvalPOSTB(:,idx_inhib)', subject(subj).zvalPOSTALL(:,idx_inhib)'};
%         tdcs(tdcs_count).psd_exc = {subject(subj).zvalSHAM(:,idx_exc)', subject(subj).zvalSTIMA(:,idx_exc)', subject(subj).zvalSTIMB(:,idx_exc)', ...
%             subject(subj).zvalPOSTA(:,idx_exc)', subject(subj).zvalPOSTB(:,idx_exc)', subject(subj).zvalPOSTALL(:,idx_exc)'};
%         tdcs(tdcs_count).psd_na = {subject(subj).zvalSHAM(:,idx_tdcsna)', subject(subj).zvalSTIMA(:,idx_tdcsna)', subject(subj).zvalSTIMB(:,idx_tdcsna)', ...
%             subject(subj).zvalPOSTA(:,idx_tdcsna)', subject(subj).zvalPOSTB(:,idx_tdcsna)', subject(subj).zvalPOSTALL(:,idx_tdcsna)'};
%         
%         % build a 4D matrix with rows= mean of psd by subject across selected channels (n), cols=freqs, z=inhib or exc
%         % or na and 4th-d=section (6 total)
%         for section = 1:length(label)
%             mean_tdcs(tdcs_count, :, 1, section) = mean(tdcs(tdcs_count).psd_inhib{1,section}); 
%             mean_tdcs(tdcs_count, :, 2, section) = mean(tdcs(tdcs_count).psd_exc{1,section});
%             mean_tdcs(tdcs_count, :, 3, section) = mean(tdcs(tdcs_count).psd_na{1,section});
% 
%             tmp_inhib(:,:,section)  = tdcs(tdcs_count).psd_inhib{1,section};
%             tmp_exc(:,:,section)  = tdcs(tdcs_count).psd_exc{1,section};
%             tmp_na(:,:,section)  = tdcs(tdcs_count).psd_na{1,section};
%         end
% 
%         tdcs_stats_inhib = [tdcs_stats_inhib; tmp_inhib];
%         tdcs_stats_exc = [tdcs_stats_exc; tmp_exc];
%         tdcs_stats_na = [tdcs_stats_na; tmp_na];
%        
%         case "tACS"
%         tacs_count = tacs_count+1;
% 
%         tacs_targ = tmp(abs(str2double(tmp.simEF)) > thrEF, ["channel", "simEF"]);
%         [~, idx] = sort(abs(str2double(tacs_targ.simEF)), 'descend');
%         tacs_targ_chan = table2array(tacs_targ(idx(1:n), "channel")); clearvars idx
%         idx_targ = ismember(subject(subj).chan, tacs_targ_chan);
%     
%         tacs_na = tmp(abs(str2double(tmp.simEF)) < thrEF, ["channel", "simEF"]);
%         [~, idx] = sort(abs(str2double(tacs_na.simEF)), 'ascend');
%         tacs_na_chan = table2array(tacs_na(idx(1:n), "channel")); clearvars idx
%         idx_tana = ismember(subject(subj).chan, tacs_na_chan);