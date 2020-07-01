function [T_chnsarea] = chanInf(chantbl_dbs, nGM, file_GMChnsarea, pdcondition)
% extract channel inf table
%
%   Output:
%       T_chnsarea: table of channel inf


nDBS = height(chantbl_dbs);

% initial the attitudes of chanarea table T_chnsarea: 
% chni_vec, electype, brainareas, notes, recordingchn
chni_vec = uint8([1: nGM+nDBS]);
electype = cell(nGM+nDBS,1);
brainareas = cell(nGM+nDBS,1);
notes = cell(nGM+nDBS,1);
recordingchn = zeros(nGM+nDBS,1);


% electype
electype(1:nGM,1) = {'Gray Matter'};
electype(nGM+1: nGM+nDBS,1) = {'DBS'};



% deal with Gray Matter
brainareas(1:nGM,1) = {''};
T = readtable(file_GMChnsarea);

if strcmp(pdcondition, 'normal')
    T.channels = T.channels_normal;
end

if strcmp(pdcondition, 'mild')
    T.channels = T.channels_mild;
end

for i = 1 : length(T.brainarea)
    area = T.brainarea{i};
    tmpcell = split(T.channels{i}, ',');
    
    for j = 1 : length(tmpcell)
        chn = str2num(char(tmpcell{j}));
        brainareas(chn,1) = {area};
        
        clear chn
    end
end
recordingchn(1:nGM) = [1:nGM];
notes(1:nGM,1) = {''};

% deal with the DBS channel table
brainareas(nGM+1 : nGM+nDBS,1) = chantbl_dbs.area;
notes(nGM+1: nGM+nDBS,1) = chantbl_dbs.elecchn;
recordingchn(nGM+1:nGM+nDBS) = [1:nDBS];

% channel information table
T_chnsarea = table;
T_chnsarea.chni = chni_vec';
T_chnsarea.brainarea = brainareas;
T_chnsarea.recordingchn = recordingchn;
T_chnsarea.electype = electype;
T_chnsarea.notes = notes;
end