function unwanted_DBS = unwanted_DBS_extract(animal)
if strcmpi(animal,'jo')
    unwanted_DBS = {'stn4-5', 'stn5-6', 'stn6-7', 'gp0-1'};
end
if strcmpi(animal,'kitty')
    unwanted_DBS = {'stn3-4', 'stn4-5', 'stn5-6', 'stn6-7', 'gp6-7'};
end
if strcmpi(animal,'pinky')
    unwanted_DBS = {'stn0-1', 'stn1-2', 'stn2-3', 'stn6-7', 'gp0-1'};
end
if strcmpi(animal,'bug')
    unwanted_DBS = {'stn5-6', 'stn6-7'};
end
end