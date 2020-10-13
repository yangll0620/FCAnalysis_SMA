function T_chnsarea = combine_GMDBSChns(T_GMchnsarea_Depth, T_DBSchnsarea)
% 	combine T_chnsarea_GM_Depth and T_chnsarea_DBS (have depth or no depth variable)
%
%

if ~any(strcmp('depth', T_DBSchnsarea.Properties.VariableNames))
    % add depth variable if T_chnsarea_DBS doesn't have
    
    T_DBSchnsarea = addvars(T_DBSchnsarea, NaN(height(T_DBSchnsarea), 1), 'NewVariableNames', 'depth', 'After', 'recordingchn');
end

T_chnsarea = vertcat(T_GMchnsarea_Depth,T_DBSchnsarea);

% adjust chi
T_chnsarea.chni = [1:height(T_chnsarea)]';

end