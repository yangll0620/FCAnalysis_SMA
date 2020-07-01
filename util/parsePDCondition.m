function condition = parsePDCondition(dateofexp, animal)
% parsePDCondtion parse the PD condition
%   parsePDCondtion_Pinky() return the PD condition (string) of one
%   particular experimental date
% 
% Input:
%   dateofexp: the experimental date in datenum format (e.g datenum(datetime('2017-07-18')))
% 
% Output:
%   condition (string): the PD condition on the dateofexp 
%                       one of {'normal', 'tomild', 'mild', 'tomoderate', 'moderate'}
%


% assign the norma, mild and moderate period
animal = lower(animal);
if strcmp(animal, 'pinky')
    period_normal = [datenum(datetime('2017-01-13')) datenum(datetime('2017-04-28'))];
    period_mild = [datenum(datetime('2017-9-12')) datenum(datetime('2017-11-27'))];
    period_moderate = [datenum(datetime('2018-5-18')) datenum(datetime('2018-12-05'))];
end
if strcmp(animal, 'bug')
    period_normal = [datenum(datetime('2017-07-18')) datenum(datetime('2018-12-23'))];
    period_mild = [datenum(datetime('2019-1-25')) datenum(datetime('2019-9-5'))];
    period_moderate = [datenum(datetime('2020-5-26'))];
end



% identify the NHP PD condition
condition = '';
if dateofexp >= period_normal(1) && dateofexp <= period_normal(2)
    condition  = 'normal';
end

if dateofexp > period_normal(2) && dateofexp < period_mild(1)
    condition = 'tomild';
end

if dateofexp >= period_mild(1) && dateofexp <= period_mild(2)
    condition = 'mild';
end

if dateofexp > period_mild(2) && dateofexp < period_moderate(1)
   condition = 'tomoderate';
end 

if dateofexp >= period_moderate(1)
    condition = 'moderate';
end


