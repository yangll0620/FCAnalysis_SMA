function notAOI_chns = notInterested_chns_extract(animal)
if strcmpi(animal,'jo')
    notAOI_chns = {};
end
if strcmpi(animal,'kitty')
    notAOI_chns = {};
end
if strcmpi(animal,'pinky')
    notAOI_chns = {'lCd', 'lSMA', 'lVA', 'lVLo', 'lVPLo', 'rMC', 'rSMA', 'rVA', 'rVLo', 'rVPLo'};
end
if strcmpi(animal,'bug')
    notAOI_chns = {};
end
end