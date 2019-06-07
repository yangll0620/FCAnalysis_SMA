function conv_tbl2struct(animal)
% python could not read the matlab table objects, thus convert the
% idxtbl_event to idxstruct_event and idx_columns

if nargin < 1
    animal = 'Pinky';
end

if isunix
    googledrive = fullfile('/home', 'lingling', 'yang7003@umn.edu');
end

if ispc
    googledrive = fullfile('F:', 'yang7003@umn');
end
folder = fullfile(googledrive, 'NMRC_umn','Projects','FCAnalysis','data', animal);
loadfolder = fullfile(folder,'epochs_reorg_filtered');
savefolder = loadfolder;
files = extractfield(dir(loadfolder),'name');
for i = 4: length(files)
    file = files{i};
    load(fullfile(loadfolder, file), 'idxtbl_event');
    idxevent = idxtbl_event{:,:};
    idx_varNames = idxtbl_event.Properties.VariableNames;
    save(fullfile(savefolder, file), 'idx_varNames','idxevent','-append');
end