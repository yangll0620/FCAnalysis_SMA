function cond_cell = cond_cell_extract(animal)
if strcmpi(animal, 'bug')
    cond_cell = {'normal', 'mild'};
end
if strcmpi(animal, 'jo')    
    cond_cell = {'normal', 'mild', 'moderate'};
end

if strcmpi(animal, 'kitty') % Kitty not have mild
    cond_cell = {'normal', 'moderate'};
end

if strcmpi(animal, 'pinky')
    cond_cell = {'normal', 'mild', 'moderate'};
end