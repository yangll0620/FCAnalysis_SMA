function [ vec_segsIndex ] = get_segIndex( states,  min_samples )
% get_segIndex 
%       get segment indices based on state and min_samples. Only get the
%       segment with state == 1 and duration is longer or equall to
%       min_samples
%                  
%       for example: state_in = [1 1 1 1 0 0 0 1 1], min_samples = 2
%                    =>   vec_segsIndex = [1 4; 8 9]
%
%       Args:
%           states: 1 * ntemp vector
%           
%           min_samples: an integer scalar
%
%       Return:
%           vec_segsIndex: nSegs * 2 (idx_start, idx_end)



%%


% find the first index at where state == 1
i = 1;
while(i <=length(states) &&states(i) == 0 )
    i =i+1;
end
idx_state1Start = i;

% extract all the seg start and end indices where state ==1 and duration >= min_samples
vec_segsIndex = [];
for i= idx_state1Start + 1: length(states)
    
    % transition from state 1 to state 0
    if states(i) == 0 && states(i-1) ==1
        idx_state1End = i -1;
        
        % the state 1 segment duration > min_samples
        if idx_state1End - idx_state1Start >= min_samples
            vec_segsIndex = cat(1, vec_segsIndex, [idx_state1Start, idx_state1End]);
        end
    end
    
    % transition from state 0 to state 1
    if states(i) == 1 && states(i-1) == 0
        idx_state1Start = i;
    end
    
    % the last state an its state is 1
    if i == length(states) && states(i) == 1
        idx_state1End = i;
        
        % the state 1 segment duration > min_samples
        if idx_state1End - idx_state1Start >= min_samples
            vec_segsIndex = cat(1, vec_segsIndex, [idx_state1Start, idx_state1End]);
        end
    end
end

