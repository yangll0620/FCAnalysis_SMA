function T_chnsarea = addDepths2chnsarea(T_chnsarea, T_chanDepth, dateofexp)
% add a new variable depths to the table T_chnsarea based on the information T_chanDepth 
% 
%   Arg:
%      dateofexp (datetime): the date of exp, e.g. dateofexp = datetime('022719', 'InputFormat', 'MMddyy'); 
%       
%      T_chnsarea: the original T_chnsarea 
%                  output of function chanInf() 
%
% 
%      T_chanDepth: table for channel depth in each dateofexp with brain area as varnames, 
%                   output of T_chanDepth_extract()
%       E.g.
%           T_chanDepth(1:3,1:5)
% 
%               3×5 table
% 
%                   Var1       M1      M1_1     M1_2      M1_3 
%                 ________    _____    ____    ______    ______
% 
%                      NaT       63     55         73        70
%                 02/27/19    26.25     50     33.125    20.625
%                 03/05/19    26.25     50     33.125    20.625

%
%  Return:
%       T_chnsarea: T_chnsarea with added variable depths
%                   [] if there is no row or more than two rows for the exp date

%% Code Start Here

% extract the row index for the exp date
idx_row = find(T_chanDepth{:, 1} == dateofexp);

% there is no row or more than two rows for the exp date
if isempty(idx_row) || length(idx_row)> 1
    if isempty(idx_row)
        disp([datestr(dateofexp) ': idx_row is empty']);
    end
    
    if length(idx_row)>1
        disp([datestr(dateofexp) ': idx_row has more than 1 row!' ]);
    end
    
    T_chnsarea = [];
    
    return;
end



% add depth 
% extract T_chanDepth__dateofexp: the channel numbers (first row) and depth (second row) for the exp date
T_chanDepth_dateofexp = [T_chanDepth(1, :); T_chanDepth(idx_row, :)];
T_chnsarea.depths = nan(height(T_chnsarea),1);
for widthi = 2: width(T_chanDepth_dateofexp) % ignore the first date column
    
    recordingchn = T_chanDepth_dateofexp{1, widthi};
    
    % find the row index for recordingchn
    idx = find(T_chnsarea.recordingchn == recordingchn);
    
    % check the brain area at recordingchn from T_chanDepth_dateofexp and T_chansarea are the same
    areastr_Tdepth = T_chanDepth_dateofexp.Properties.VariableNames{widthi}; % areastr_Tdepth = 'M1' or 'M1_1', 'M1_2'
    area_Tchnsarea = T_chnsarea.brainarea{idx};
    
    if ~contains(areastr_Tdepth, area_Tchnsarea(find(~isspace(area_Tchnsarea))))
        disp([datestr(dateofexp) ', recordingchn = ' num2str(recordingchn) ': brainarea for T_chnsarea = ' area_Tchnsarea ' , T_chanDepth = ' areastr_Tdepth])
        
        newstr = split(areastr_Tdepth, '_');
        
        T_chnsarea.brainarea{idx} = newstr{1};

        disp(['change T_chanares.brainarea to ' newstr{1}])
    end
    
    T_chnsarea.depths(recordingchn) = T_chanDepth_dateofexp{2, widthi};
    
    clear recordingchn
    
end
end