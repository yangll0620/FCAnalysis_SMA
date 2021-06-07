function func_replaceAllNHPsFunc()
str = input('input function name (without .m): ','s'); 
filename = [str '.m'];

% Code starting here
folder_current = pwd;

% animal
if ismac
    % Code to run on Mac platform
elseif isunix
    % Code to run on Linux platform
    
    [fi, j] = regexp(folder_current, ['NHPs', '/', '[A-Za-z]*']);
elseif ispc
    % Code to run on Windows platform
    
    [fi, j] = regexp(folder_current, ['NHPs', '\\', '[A-Za-z]*']);
else
    disp('Platform not supported')
end
animal = folder_current(fi + length('NHPs') + 1:j);


NHPs = {'Jo', 'Kitty', 'Bug', 'Pinky'};
disp(['coping '  fullfile(folder_current, filename)])
for ni = 1 : length(NHPs)
    if strcmpi(NHPs{ni}, animal)
        continue;
    end
    folder_des = strrep(folder_current, animal, NHPs{ni});
    status = copyfile(fullfile(folder_current, filename), fullfile(folder_des, filename));
    if status
        disp(['copied ' filename ' into ' folder_des])
    end
    clear folder_des status
end


