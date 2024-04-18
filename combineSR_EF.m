
clear
spikerate = readtable("C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\PS2_OFFICIAL\ChanInfo_PS2.xlsx", "sheet", "SR");

spikerate.subject = string(spikerate.subject);
idx = strcmp(spikerate.subject,"");
spikerate(idx,:) = [];

[subjs,ia,ic] = unique(spikerate.subject, 'stable');
ia(end+1) = length(spikerate.subject)+1;
%spikerate.channel = string(spikerate{:,2});


dir_data = 'C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani';
mono = 0;
EFtable = EFtable(dir_data, 'monopolar', mono);

for subj = 1:length(ia)-1
    clearvars subchan subtable i simEF
    subchan = spikerate.Channel(ia(subj):ia(subj+1)-1,:);
    subchan = string(subchan);
    subtable = EFtable{strcmp(EFtable(:,1), subjs(subj)),2}; %channel list of EFtable corresponding to current subject
    for i = 1:length(subchan)
        simEF = subtable{strcmp(subtable.Contact, subchan(i)),8};
        if ~isempty(simEF)
            spikerate.simEF(i+ia(subj)-1) = simEF;
        else spikerate.simEF(i+ia(subj)-1) = NaN;
        end
    end
    
end

        