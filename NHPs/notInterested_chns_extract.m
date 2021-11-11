function noisy_chns = notInterested_chns_extract(animal)
if strcmpi(animal,'jo')
    noisy_chns = {};
end
if strcmpi(animal,'kitty')
    noisy_chns = {};
end
if strcmpi(animal,'pinky')
    noisy_chns = {'lCd', 'lSMA', 'lVA', 'lVLo', 'lVPLo', 'rMC', 'rSMA', 'rVA', 'rVLo', 'rVPLo'};
end
if strcmpi(animal,'bug')
    noisy_chns = {};
end
end