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
            aux = h5read(h5name,['/' info.Datasets(j).Name]);
            DataE(fi,:,:) = aux;
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
for ic = 1 : size(Data,1)
    100*ic/size(Data,1)
    capteur(ic).Name = capName{ic};
    capteur(ic).Pos = Pos(ic,:);
    aux = squeeze(Data(ic,:,:));
    capteur(ic).Time = aux(1,1:subs:end);
    aux(1,:) = [];

    if any(labels == "EnergyP    1")
        capteur(ic).EnergyP = aux(1,1:subs:end);
        aux(1,:) = [];
    end
    if any(labels == "EnergyK    1")
        capteur(ic).EnergyK = aux(1,1:subs:end);
        aux(1,:) = [];
    end
    if any(labels == "Eps Vol    1")
        capteur(ic).EpsVol = aux(1,1:subs:end);
        aux(1,:) = [];
    end
    if any(labels == "Displ      1")
        capteur(ic).Displ = aux(1:3,1:subs:end);
        aux(1:3,:) = [];
    end
    if any(labels == "Veloc      1")
        capteur(ic).Veloc = aux(1:3,1:subs:end);
        aux(1:3,:) = [];
    end
    if any(labels == "Accel      1")
        capteur(ic).Accel = aux(1:3,1:subs:end);
        aux(1:3,:) = [];
    end
    if any(labels == "Pressure   1")
        capteur(ic).Pressure = aux(1,1:subs:end);
        aux(1,:) = [];
    end  
    if any(labels == "Eps Dev    1")
        capteur(ic).EpsDev = aux(1:6,1:subs:end);
        aux(1:6,:) = [];
    end
    if any(labels == "Stress Dev 1")
        capteur(ic).StressDev = aux(1:6,1:subs:end);
        aux(1:6,:) = [];
    end
    if any(labels == "Eps Dev Pl 1")
        capteur(ic).EpsDevPl = aux(1:6,1:subs:end);
        aux(1:6,:) = [];
    end
    if any(labels == "DUDX       1")
        capteur(ic).DUDX = aux(1:9,1:subs:end);
        aux(1:9,:) = [];
    end
    if any(labels == "GradLambda 1")
        capteur(ic).GradLambda = aux(1:3,1:subs:end);
        aux(1:3,:) = [];
    end
    if any(labels == "GradMu     1")
        capteur(ic).GradMu = aux(1:3,1:subs:end);
        aux(1:3,:) = [];
    end
    if any(labels == "EnergyD    1")
        capteur(ic).EnergyD = aux(1:3,1:subs:end);
        aux(1:3,:) = [];
    end
end

save([folder filesep matname '.mat'],'capteur','DataE','-v7.3')