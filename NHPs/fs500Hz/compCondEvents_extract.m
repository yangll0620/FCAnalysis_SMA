function [tbl_compConds, tbl_compEvents]= compCondEvents_extract(animal, varargin)

p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
parse(p,varargin{:});
codesavefolder = p.Results.codesavefolder;

% copy code to savefolder if not empty
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end

% extract tbl_compConds
cond_cell = cond_cell_extract(animal);
nconds = length(cond_cell);
tbl_compConds = table;
tblVarNames_cond = {'baseCond', 'compCond'};
for baseci = 1 : nconds -1
    cond_base = cond_cell{baseci};
    for compci = baseci+1 : nconds
        cond_comp = cond_cell{compci};
        tbl_compConds = [tbl_compConds; cell2table({cond_base, cond_comp}, 'VariableNames', tblVarNames_cond)];
        clear cond_comp
    end
    clear cond_base
end

% extract tbl_compConds
EventPhases = SKT_eventPhases_extract(animal);
nevents = length(EventPhases);
tbl_compEvents = table;
tblVarNames_event = {'baseEvent', 'compEvent'};
for baseei = 1 : nevents -1
    event_base = EventPhases{baseei};
    for compei = baseei+1 : nevents
        event_comp = EventPhases{compei};
        tbl_compEvents = [tbl_compEvents; cell2table({event_base, event_comp}, 'VariableNames', tblVarNames_event)];
        clear event_comp
    end
    clear event_base
end

end