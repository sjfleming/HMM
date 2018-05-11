function obs_prob = emissions_prob(measurements, states, p_noise)
% emissions_probs(measurement, level_mean, level_stdv, p_noise)
% returns the probability of a bunch of measured currents given the 
% expectation.
% uses the KL-divergence to calculate distance
% (https://en.wikipedia.org/wiki/Kullback-Leibler_divergence)
% Inputs:
% 'measurements': a cell array of current arrays for different states
% 'states' is an M-element cell array that contains three fields:
%       'level_mean' level's mean current
%       'level_stdv' standard deviation of level current
%       'stdv_mean' standard deviation of level's mean current
% 'p_noise' the a priori probability of observing a junk noise value
% Returns:
% 'obs_prob': a matrix of states by measurements with their probabilities
% Stephen Fleming
% 4/5/18
    
    % finds the overlap between current measurements given an 
    % expected level mean and standard deviation.
    
    obs_prob = nan(numel(states),numel(measurements));
    
    for s = 1:numel(measurements)
        n = numel(measurements{s});
        p = ones(n,1) / max(1,n); % data probs
        q = cellfun(@(x) normpdf(measurements{s},x.level_mean,x.level_stdv), states); % model probs
        obs_prob(s,:) = -1*sum(p.*log(q./p)) + p_noise;
    end
        
    % normalize if needed
    m = max(max(obs_prob));
    if m>1
        obs_prob = obs_prob / m;
    end
    
end