function condition = parsePDCondition_Pinky(dateofexp)
% parsePDCondtion_Pinky parse the PD condition
%   parsePDCondtion_Pinky() return the PD condition (string) of one
%   particular experimental date
% 
% Input:
%   dateofexp: the experimental date in datenum format
% 
% Output:
%   condition (string): the PD condition on the dateofexp 
%                       one of {'normal', 'tomild', 'mild', 'tomoderate', 'moderate'}
%

period_normal = [datenum(datetime('2017-01-13')) datenum(datetime('2017-04-28'))];
period_mild = [datenum(datetime('2017-9-12')) datenum(datetime('2017-11-27'))];
period_moderate = [datenum(datetime('2018-5-18')) datenum(datetime('2018-12-05'))];

% identify the NHP PD condition
if dateofexp >= period_moderate(1)
    condition = 'moderate';
else
    if dateofexp >= period_mild(2)
        condition = 'tomoderate';
    else 
        if dateofexp >= period_mild(1)
            condition = 'mild';
        else
            if dateofexp >= period_normal(2)
                condition = 'tomild';
            else
                if dateofexp >= period_normal(1)
                    condition  = 'normal';
                end
            end
        end
    end
end