%% rename files in a folder by keeping their subject code

% delphos
folder = "\\dynaserv\Galvani_ps2\montages\delphos";
myfiles  = dir(strcat(folder, "\*.mtg"));

for i=1:length(myfiles)
    subj = extractBefore(myfiles(i).name, '_');
    newname = strcat(subj, "_bipolar-raw.mtg");
    movefile(strcat(folder, '\', myfiles(i).name), strcat(folder, '\', newname))
end

% FC
folder = "\\dynaserv\Galvani_ps2\montages\FC";
myfiles  = dir(strcat(folder, "\*.mtg"));

for i=1:length(myfiles)
    subj = extractBefore(myfiles(i).name, '_');
    newname = strcat(subj, "_bipolar-selection.mtg");
    movefile(strcat(folder, '\', myfiles(i).name), strcat(folder, '\', newname))
end

% Delphos results 
folder = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\Delphos\raw\**\*.mat";
output_fold = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\Delphos\all-together";
myfiles  = dir(folder);
tochange = contains(string({myfiles.folder}), 'sub', 'IgnoreCase', true);

for i=1:length(myfiles)
    subj = extractAfter(myfiles(i).folder, '_');
    if  ~contains(myfiles(i).name, 'sub')
        newname = strcat(subj, '_', myfiles(i).name);
        copyfile(strcat(myfiles(i).folder, '\', myfiles(i).name), strcat(output_fold, '\', newname))
    else 
        copyfile(strcat(myfiles(i).folder, '\', myfiles(i).name), strcat(output_fold, '\', myfiles(i).name))
    end
end


