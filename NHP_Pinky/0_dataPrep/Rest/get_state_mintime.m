function [state_vec] = get_state_mintime( state_in, min_samples )
%GET_STATE_MINTIME 
%       get min state vector based on state_in and min_samples.
%       state_vec is set to be 1/0 if state_in is 1 and the state 1/0 last longer than min_samples, others -1
%                  
%       for example: state_in = [1 1 1 1 0 0 0 1 1], min_samples = 3
%                    =>   state_vec = [1 1 1 1 0 0 0 -1 -1]
%
%       Args:
%           state_in: 1 * ntemp vector
%           
%           min_samples: an integer scalar
%
%       Return:
%           state_vec: 1 * ntemp vector


%%
state_vec = -1*ones(size(state_in));
idx_start = 1; % start index of one state segment
for k = 2:length(state_in)
    
    % state change to a different state (change point) or the last state 
    if ( state_in(k) ~=  state_in(k-1) || k==length(state_in) )
        
        % end index of the last state segment
        idx_end = k-1;
        
        % the last state segement longer than min_elapsed_samples
        if idx_end-idx_start +1 >= min_samples
            
            % previous state is 1, state_vec for this state is assigned 1
            if state_in(idx_start) ==1
                state_vec(idx_start:k-1) = ones(1, k - idx_start )  ;
            end
            
            % previous state is 0, state_vec for this state is assigned 0
            if state_in(idx_start) ==0
                state_vec(idx_start:k-1) = zeros(1, k - idx_start );
            end
        end
        
        idx_start = k;
    end
end

end

