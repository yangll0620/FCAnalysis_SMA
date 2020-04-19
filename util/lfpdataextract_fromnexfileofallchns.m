function [lfpdata, fs] = lfpdataextract_fromnexfileofallchns(nexfolder, file_prefix, chns)
%% extract lfp data from nex files
%   Inputs:
%          nexfolder: the folder storing all the nex files 
%               (e.g /home/lingling/root2/Animals2/Pinky/Recording/Processed/DataDatabase/Pinky_011918/LFP/Block-2)
%          file_prefix the prefix of the nex files 
%               (e.g Pinky_GrayMatter_eyetracking_DT1_011918_Block-2_LFPch)
%          chns: all the channel numbers to be extracted
 
for i = 1: length(chns)
    
    filename = [file_prefix num2str(chns(i)) '.nex'];
    
    % read nex file data
    [nexdata] = readNexFile(fullfile(nexfolder, filename));
    
    % extract the number of the structure containing LFPchn* data
    name_list = extractfield(cell2mat(nexdata.contvars),'name');
    
    % nexlfp.contvars name == 'LFPch1', or 'MUAch1'
    i_lfp = find(contains(name_list, 'LFP'));
    
    % lfp data
    data = nexdata.contvars{i_lfp}.data;
    
    if i == 1 
        % sample rate
        fs = nexdata.contvars{i_lfp}.ADFrequency;
        
        % lfpdata is for the lfp data of all the chns
        lfpdata = zeros(length(data), length(chns));
   
    else
        % samping frequency is different
        if fs ~= nexdata.contvars{i_lfp}.ADFrequency
            chni = chns(i);
            disp([nexfolder 'fs = ' num2str(fs) ',sampling frequency is different for chni = ' num2str(chni)]);
            lfpdata = [];
            fs = [];
            break;
        end
    end
    
    % assign to lfpdata
    lfpdata(:, i) = data;
    
    clear filename nexdata name_list i_lfp data i
end