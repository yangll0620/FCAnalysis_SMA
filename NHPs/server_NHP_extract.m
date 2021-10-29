function server_NHP = server_NHP_extract(animal)

if strcmpi(animal, 'Jo') || strcmpi(animal, 'Kitty')
    server_NHP = fullfile('Z:', 'root2', 'Animals');
end