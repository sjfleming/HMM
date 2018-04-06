function obs_prob = emission_prob(measurement, states, p_noise)
% emission_probs(measurement, level_mean, level_stdv, p_noise)
% returns the 2-tailed p-value for a measured current given the expectation
% Inputs:
% 'measurement': a column vector of current values
% 'measured_current_range' row 2-vector of of the entire current range
% 'states' is an M-element cell array that contains three fields:
%       'level_mean' level's mean current
%       'level_stdv' standard deviation of level current
%       'stdv_mean' standard deviation of level's mean current
% 'p_noise' the a priori probability of observing a junk noise value
% Stephen Fleming
% 4/5/18
    
    % finds the two-tailed p-value for a current measurement given an 
    % expected level mean and standard deviation.
    
    if size(measurement,2)~=1 % flip row vector to column vector if needed
        measurement = measurement';
    end
    
    % (see doc erf, look at the CDF of a Gaussian, multiply by 2)
    % adding in p_noise as if it integrates to p_noise over (-Inf, Inf)
    obs_prob = cell2mat(cellfun(@(s) 1 + erf(-0.7071*abs(measurement - s.level_mean)/s.level_stdv), ...
        states, 'uniformoutput', false))' + p_noise;
    
    % normalize if needed
    m = max(max(obs_prob));
    if m>1
        obs_prob = obs_prob / m;
    end
    
end