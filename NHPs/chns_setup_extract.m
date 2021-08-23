function [nM1, nGM, nSTN, nGP] = chns_setup_extract(animal)

if strcmpi(animal, 'Pinky')
    nM1 = 96; nGM = 32; nSTN = 8; nGP = 8;
end

if strcmpi(animal, 'Jo')
    nM1 = 96; nGM = []; nSTN = 8; nGP = 8;
end

