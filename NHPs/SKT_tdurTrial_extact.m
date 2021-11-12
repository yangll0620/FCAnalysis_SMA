function [tdur_trial_normal, tdur_trial_mild, tdur_trial_moderate] = SKT_tdurTrial_extact(animal)

tdur_trial_normal = [];
tdur_trial_mild = [];
tdur_trial_moderate = [];

if strcmpi(animal, 'bug')
    tdur_trial_normal = [-1 1];
    tdur_trial_mild = [-1 1];
end
if strcmpi(animal, 'jo')
    tdur_trial_normal = [-1 1];
    tdur_trial_mild = [-1 1];
    tdur_trial_moderate = [-1 1];
    
end

if strcmpi(animal, 'kitty') % Kitty not have mild
    tdur_trial_normal = [-1 1];
    tdur_trial_moderate = [-1 1];
end

if strcmpi(animal, 'pinky')
    tdur_trial_normal = [-1 1];
    tdur_trial_mild = [-1 1];
end