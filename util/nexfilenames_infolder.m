function [file_prefix, chns] = nexfilenames_infolder(nexfolder)
%% return the prefix of file and chns


%%
filenames = extractfield(dir(nexfolder), 'name');

% match channel file
match = cellfun(@(x) ~isempty(regexp(x, ['LFPch[0-9]*.nex'],'match')),filenames,'UniformOutput',false); 
match = cell2mat(match);

nexfilenames = filenames(match);

file_prefix = '';
chns = [];

% parse and sort the number of channel
for filei = 1 : length(nexfilenames)
    nexfilename = nexfilenames{filei}; 
    tmp = char(regexp(nexfilename, 'ch[0-9]*+.nex', 'match'));
    chns(filei) = str2num(tmp(length('ch')+1: end-length('.nex')));
    
    % get the prefix of the each nex file name
    tmp = split(nexfilename, [num2str(chns(filei)) '.nex']);
    
    if filei == 1
        file_prefix = tmp{1};
    else
        if ~strcmp(file_prefix, tmp{1})
            file_prefix = '';
            disp('file_prefix is not consistent')
            break;
        end
    end
end
chns = sort(chns);