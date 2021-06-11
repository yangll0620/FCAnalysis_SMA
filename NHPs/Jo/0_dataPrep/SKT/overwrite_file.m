function output = overwrite_file(fileinfo)

% Construct a questdlg with three options
choice = questdlg(['Would you like to overwrite existing file? ' fileinfo], ...
	'Overwrite File?', ...
	'Yes','No','No');
% Handle response
switch choice
    case 'Yes'
        output = 1;
    case 'No'
        output = 0;

end

