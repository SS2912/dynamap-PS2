function EFtable_melt = EFtable(dir_data, varargin)

% function to obtain list of bipolar and monopolar channels with simulated
% EF from Neuroelectrics in the context of PS2. One file for patient, will
% be put in a cell with patient code (first part of file name, before "_")

% dir_data: directory of the EF tables (general directory with all the subfolders containing the tables)
% optional: 'monopolar': do you want to have one EF value per channel (monopolar =1) or do you want the bipolar channels with mean EF for that couple (monopolar=0, default)?

% Debug and example:
% dir_data = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani";
% mono = 1;
format long g
subdirs = strcat(dir_data, '\**\*.txt');
files = dir(subdirs);
stringToBeFound = 'Efield_at_seeg_contacts';
EFtable = {};
mono = 0; %default: bipolar list
for ii = 1:2:nargin-1
        if strcmp('monopolar', varargin{ii})
            mono = varargin{ii+1}; 
        end
end
    
idx=0;

for i=1:length(files)
    if contains(files(i).name, stringToBeFound)
        idx=idx+1;
        EFtable{idx,1} = extractBefore(files(i).name,"_");
        cd(files(i).folder);
        temp = readtable(files(i).name);
        
        temp.Contact = string(temp.Contact);
        monopolar = strcat(temp.Electrode, temp.Contact);
        bipolar   = strcat(monopolar(1:end-1), '-', monopolar(2:end));   %will also create elements of different electrodes (ex:"OF14-CC1"), but they'll be deleted when comparing this list to the montage
        if mono
            temp(:,2) = table(monopolar);
            EFtable{idx,2} = temp(:,2:end);
        else
            temp(1:end-1,2) = table(bipolar);
            for j=1:size(temp,1)-1
                temp(j,9) = table(mean(temp.E_field(j:j+1)));
            end
           
            EFtable{idx,2} = temp(1:end-1,2:end);
        end
        
        clearvars monopolar bipolar temp
    end
    
end

%% create a big table with all the channels and subject code for each channel

row = 1;
clearvars EFtable_melt
EFtable_melt = {};
for i = 1:length(EFtable)
    l = size(EFtable{i,2}, 1);
    EFtable_melt(row:row+l-1,1) = repelem(EFtable(i,1), l)';
    EFtable_melt(row:row+l-1,2:9) = table2cell(EFtable{i,2});
    row = row+l;
end


EFtable_melt = cell2table(EFtable_melt, 'VariableNames',{'subject','channel', 'x', 'y', 'z', 'isgrey', 'matter', 'brain_area', 'simEF'});

end
        
       

