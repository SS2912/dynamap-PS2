%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                SPIKE RATE ANALYSIS - PS2                   %%%%%%%
%%%%%%%                  S. Simula - Feb 2023                      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% functions needed: 
% - dyn_read_mtg.m,
% - mean_rate_delphos.m, 


%% 1a. prepare to read the .mat files in the folders 
% to change with your directory if different than this one. 
clear

% input directory (with delphos results): 
input_folder  = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\Delphos\mat_threshold80";
sections = {"baseline", "sham", "postA", "postB", "postall"};

detections = struct();

%% 1b. read the montages of each patient, in which you will insert the calculated SR 
% (step which is necessary here cause delphos only keeps the non-zero detections channels in its results tables)

SR_diff_table = table();
mtgs = dir('\\dynaserv\Galvani_ps2\montages\delphos\*.mtg');

results_all = struct();
delphos_resultsALL = table();

for i = 1:size(mtgs,1) % for each subject
    delphos_subtab = table();
    sub = string(extractBetween(mtgs(i).name, '-', '_'));
    montage_struct = dyn_read_mtg(strcat(mtgs(i).folder, '\', mtgs(i).name));
    results_all(i).subj = sub(1);

    % create table with bipolar channels and sub code where you'll insert SR info
    for j = 1:length(montage_struct)
        delphos_subtab.subj(j)= sub(1);  % sub(1) because in old file names there were different string delimited by - and _ -> I only took the first one (subj code)
        delphos_subtab.chan(j) = strcat(string(montage_struct(j).name), '-', string(montage_struct(j).reference));
    end
    

    %% 2. find the same subject in delphos input folder and read + write the results onto SR_diff_table
    % 2a. for each section(folder), find the current subject and load their delphos results
    idxev = 0;

    for ii = 1:length(sections)
        curr_folder = strcat(input_folder, "\", sections{ii});
        cd(curr_folder)
        myfiles  = dir(strcat(curr_folder, "\*.mat"));
        match    = contains(string({myfiles.name}), sub(1), 'IgnoreCase', true);
       
        if sum(match)
           load(myfiles(match).name);
           windowsize    = 60;
           events_labels = string({results.markers.label});
           events_list   = cellstr(unique(events_labels));
           if size(results.labels,1)<size(results.labels,2)
              chan_list   = string(results.labels)';
           else 
              chan_list = string(results.labels);
           end
           
           % create big table where to insert matching channels from delphos results to raw channel list
           subtab_event         = table();           
           subtab_event.chan    = chan_list;
           subtab_event.section = repelem(sections{ii}, length(chan_list), 1);
           colstart             = size(subtab_event,2);

    % 2b. for each subject, and for each section, calculate spike/osc rate via the "mean_rate_delphos.m function
        for ev = 1:length(events_list)
           
           event = events_list{ev};
           detections(ev +idxev).subj                 = sub(1);
           detections(ev +idxev).section              = sections{ii};
           detections(ev +idxev).event                = event; 
           detections(ev +idxev).channels             = chan_list;  
           detections(ev +idxev).windowsize           = windowsize;
           [detections(ev +idxev).meanrate_window, ~] = mean_rate_delphos(results, windowsize, event); % mean rate per window (mean_window) and timepoints of all detections for each channel
         
           
           subtab_event(:,ev+colstart) = table(mean(detections(ev).meanrate_window, 2)*(60/windowsize));
            
%            tmp = delphos_subtab;
%           
        end  
       
        idxev = 1+idxev+length(sections);
        results_all(i).detections = detections;

        % now match the channels in delphos with the ones in raw mtg
        subtab_event.Properties.VariableNames = ["chan", "section", events_list{:}];
        for chan=1:length(delphos_subtab.chan)
            chan_match = strcmp(string(delphos_subtab.chan(chan)), string(subtab_event.chan));
            if ~isempty(subtab_event(chan_match, :))
                delphos_subtab(chan,3:10) = subtab_event(chan_match, :);
            else
                delphos_subtab(chan,3:10) = {NaN};
            end
        end

        delphos_resultsALL = [delphos_resultsALL; delphos_subtab];
        end
        
    end

    clearvars delphos_subtab detections colstart sub
end

delphos_resultsALL.Properties.VariableNames = ["subj", "chanRaw", "chanDelphos", "section", events_list{:}];


%         
%         for chan=1:size(subtab_event,1)
%         SR_subtab(strcmp(string(subtab_pre.Channels(chan)), SR_subtab.chan), 3:length(sections)+3 = subtab_pre.Spk_Rate(chan);
%         end
% subtab_event
%           for chan=1:size(subtab_pre,1)
%                 SR_subtab.SRpre(strcmp(string(subtab_pre.Channels(chan)), SR_subtab.chan)) = subtab_pre.Spk_Rate(chan);
%             end
%     % 
% 
%     
% 
% 
% 
% 
% 
% 
%             for chan=1:size(subtab_pre,1)
%                 SR_subtab.SRpre(strcmp(string(subtab_pre.Channels(chan)), SR_subtab.chan)) = subtab_pre.Spk_Rate(chan);
%             end
%         else 
%             SR_subtab.SRpre(:) = 0;
%         end
%         
%         clearvars results subtab
%         SR_subtab('VariableNames',["subj","chan", sections])
% 
%     end
% 
%     
% end
% SR_diff_table.SRdiff = SR_diff_table.SRpost - SR_diff_table.SRpre;
% 
% name = strcat('C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Thermo\spikerate','\','SRtable.csv');
% writetable(SR_diff_table,name)  

% if nsub_post < nsub_pre
%     allsubj = unique(tab_pre.sub);
%     idx = ismember(allsubj, unique(tab_post.sub), 'rows');
%     missing_sub = tab_pre(strcmp(tab_pre.sub, allsubj(~idx)),:);   
% else if nsub_post > nsub_pre
%         idx = ismember(unique(tab_post.sub), unique(tab_pre.sub), 'rows');
%         missing_sub = tab_pre(strcmp(tab_pre.sub, allsubj(~idx)),:);   
%     end
% end
% 
% %% merge the two tabs and calculate the SR difference for each channel


% In delphos github:
% mrk_chn_bln = strcmp({detection_markers_chn(:).label}, handles.markers2display{mrk_idx});
% handles.rate_markers(i,mrk_idx) = sum(mrk_chn_bln)/(handles.duration/60);
            