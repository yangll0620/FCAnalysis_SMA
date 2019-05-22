clear
varnames = {'animal', 'dateofexp', 'block', 'condition', 'therapy', 'badtrials'};
T = table({'Jo'}, {'10-13-15'}, 8, {'normal'},{''}, {[8,9,14,15,48,51,80]} ,'VariableNames',varnames);
T = [T;{'Jo','10-14-15', 10,'normal','', [1]}];

T = [T;{'Jo','11-18-15', 2,'mild','', []}];
%T = [T;{'Jo','11-18-15', 15,'mild','ldopa'}];
T = [T;{'Jo','11-19-15', 2,'mild','', []}];
%T = [T;{'Jo','11-19-15', 14,'mild','ldopa'}];
T = [T;{'Jo','11-24-15', 2,'mild','', [1]}];
%T = [T;{'Jo','11-24-15', 17,'mild','ldopa'}];

T = [T;{'Jo','01-20-16', 2,'moderate','', []}];
%T = [T;{'Jo','01-20-16', 11,'moderate','dbsgpi'}];
% T = [T;{'Jo','01-22-16', 2,'moderate','', []}]; % not successful
%T = [T;{'Jo','01-22-16', 10,'moderate','dbsgpi'}];
T = [T;{'Jo','02-04-16', 3,'moderate','', [1]}];
T = [T;{'Jo','02-16-16', 1,'moderate','', [1]}];
for i = 8: height(T) 
    animal = char(T{i, 'animal'});
    dateofexp = char(T{i, 'dateofexp'});
    block = T{i, 'block'};
    condition = char(T{i, 'condition'});
    therapy = char(T{i, 'therapy'});
    badtrials = cell2mat(T{i, 'badtrials'});
    extractlfptrial_m1dbs(animal, dateofexp , block, condition, therapy, badtrials);
    
    clear animal dateofexp block condition therapy
end

writetable(T,'lfptrialdataset.csv','Delimiter',','); 