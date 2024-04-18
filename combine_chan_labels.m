function combine_chan_labels(labelsA, varsA, labelsB, varsB, varargin)

% function to combine values and/or labels of the same channels from different
% tables. For example, combine together EI values and labels (EZ, PZ etc)
% with PSD values from a different file. The frist input table will be the
% one dispalyed first so better to have subject list there.
% input: 
% - labelsA, labelsB: full path name of excel or text files with a channel column and at least another col with labels or values
% - varsA, varsB: names of variables (columns) to be combined. options: string of column title or "all" to take all columns
% varargin: 
% - monopolar or bipolar (cuts bipolar channel list if option is monopolar)
% - sheetA, sheetB: excel sheets to read (EI tables: "Synthesis")
% - subjA, subjB: subject code if one or both files need a selection of subject
% output:
% table with combines values or labels for the channel list inputed on alphabetical order

%% read varargin

for ii = 1:2:nargin-4
        if strcmp('monopolar', varargin{ii})
            monopolar = varargin{ii+1}; 
        elseif strcmp('sheetA', varargin{ii})
            sheetA = varargin{ii+1};
        elseif strcmp('sheetB', varargin{ii})
            sheetB = varargin{ii+1};
        elseif strcmp('subjA', varargin{ii})
            subjA = varargin{ii+1};
        elseif strcmp('subjB', varargin{ii})
            subjB = varargin{ii+1};
        end
end

%% example:

% clear 
% clc
% labelsA = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\EI_tables\sub-0b12810f4f6d_TE_EIsummary.xlsx";
% varsA = ["Channel" "max_cEI" "ROI_cEI" "Structure"];
% sheetA = "Synthesis";
% 
% labelsB = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\SpikeRate_PS2.xlsx";
% varsB = "all";  %["subject" "Channel" "SR_base" "SR_sham" "SR_postA" "SR_postB" "SR_postALL" "simEF" "group"];
% sheetB = "SR";
% subjB = "0b12810f4f6d"; 

%% read tables

if exist('sheetA', 'var')
    tableA = readtable(labelsA, 'Sheet', sheetA);
else
    tableA = readtable(labelsA);
end

if exist('sheetB', 'var')
    tableB = readtable(labelsB, 'Sheet', sheetB);
else
    tableB = readtable(labelsB);
end

% change channel column to 'channel' (sometimes is Channels) and uniform
% everything to low case
idx = contains(tableA.Properties.VariableNames,'channel','IgnoreCase',true);
tableA.Properties.VariableNames(idx) = {'channel'}; clearvars idx;
idx = contains(tableA.Properties.VariableNames,'subj','IgnoreCase',true);
tableA.Properties.VariableNames(idx) = {'subject'}; clearvars idx;
%tableA.Properties.VariableNames = lower(tableA.Properties.VariableNames);

idx = contains(tableB.Properties.VariableNames,'channel','IgnoreCase',true);
tableB.Properties.VariableNames(idx) = {'channel'}; clearvars idx;
idx = contains(tableB.Properties.VariableNames,'subj','IgnoreCase',true);
tableB.Properties.VariableNames(idx) = {'subject'}; clearvars idx;
%tableB.Properties.VariableNames = lower(tableB.Properties.VariableNames);

%% select variables in tables

if ~strcmp(varsA, "all")
    varsA = lower(varsA);
    tableA = tableA(:, varsA);
end

if ~strcmp(varsB, "all") 
    varsB = lower(varsB);
    tableB = tableB(:, varsB);
end

if exist('subjA', 'var')
   tableA = tableA(tableA.subject == subjA, :);
end

if exist('subjB', 'var')
   tableB = tableB(tableB.subject == subjB, :);
end

%% combine the two tables

idx = 0;
for i = 1:size(tableA,1)
    sel = tableB(tableB.channel == string(tableA.channel(i)),:);
    if ~isempty(sel)
        idx = idx+1;
        table(idx,:) = join(tableA(i,:), sel);
    end
    clearvars sel
end

end