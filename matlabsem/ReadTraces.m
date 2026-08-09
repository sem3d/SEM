close all
clear all
clc

% load the files...
folder = '';
matname = 'case';


files = dir([folder filesep 'capteurs*.h5']);

%% load the data
subs = 1;

for ifile = 1 : numel(files)
    h5name = [files(ifile).folder filesep files(ifile).name];
    info = h5info(h5name);
    for j = 1 : numel(info.Datasets)
        j/numel(info.Datasets)
        if(strcmp(info.Datasets(j).Name, 'Variables'))
            labels = h5read(h5name, '/Variables');
        elseif(strcmp(info.Datasets(j).Name, 'Energy_Variables'))
            labelsE = h5read(h5name, '/Energy_Variables');
        elseif strcmp(info.Datasets(j).Name(end-3:end), '_pos')
            name = info.Datasets(j).Name;
            if strcmp(name, 'Energy_pos')
                %do nothing
            else
                cap = regexp(name,'_','split');
                Pos(str2double(cap{2})+1,:) = h5read(h5name,['/' name])';
            end
        elseif strcmp(info.Datasets(j).Name, 'Energy')
            fi = sscanf(files(ifile).name,'capteurs.%d.h5')+1;
            if fi == 1
                aux = h5read(h5name,['/' info.Datasets(j).Name]);
                DataE = aux;
            end
        else
            name = info.Datasets(j).Name;
            cap = regexp(name,'_','split');
            capName{str2double(cap{2})+1} = name;
            aux = h5read(h5name,['/' name]);
            Data(str2double(cap{2})+1,:,:) = aux(:,1:subs:end);
        end
    end
end
%% create capteur structure
labelsStr = strtrim(string(labels));
baseNames = regexprep(labelsStr, '\s+\d+$', '');
uniqueNames = unique(baseNames, 'stable');

for ic = 1 : size(Data,1)
    100*ic/size(Data,1)
    capteur(ic).Name = capName{ic};
    capteur(ic).Pos = Pos(ic,:);
    aux = squeeze(Data(ic,:,:));

    for iv = 1 : numel(uniqueNames)
        uName = uniqueNames(iv);
        fieldName = char(strrep(uName, " ", ""));
        idx = find(baseNames == uName);
        capteur(ic).(fieldName) = aux(idx, 1:subs:end);
    end
end

save([folder filesep matname '.mat'],'capteur','DataE','-v7.3')