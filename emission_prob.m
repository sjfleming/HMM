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
    
    % this is a 2-tailed p-value for a Cauchy distribution, which goes to
    % zero MUCH more slowly than a Gaussian
%     obs_prob = cell2mat(cellfun(@(s) 1 + 2*atan(-abs(measurement-s.level_mean)/s.level_stdv)/pi, ...
%         states, 'uniformoutput', false))' + p_noise;
    
    % maybe try chi-squared PDF for 1 degree of freedom... chiPDF( (x-u)^2/sig^2 )
    % 1/(sqrt(2)*gamma(0.5)) = 0.3989422804, k = 1
    % since x^k/2-1 is 1/sqrt(x) and x is (x-u)^2/sig^2, this is
    % sig/abs(x-u)
%     obs_prob = cell2mat(cellfun(@(s) 0.7978845608 * exp(-(measurement-s.level_mean).^2/(s.level_stdv^2)/2) ...
%         * s.level_stdv ./ abs(measurement-s.level_mean), ...
%         states, 'uniformoutput', false))' + p_noise;
    
    % try a t-distribution
    % still 2-tailed p-value from the CDF, 1 dof
%     obs_prob = cell2mat(cellfun(@(s) 2*tcdf(-abs(measurement-s.level_mean) / s.level_stdv, 1), ...
%         states, 'uniformoutput', false))' + p_noise;
    % oh, for 1 dof, this seems to be identical to the Cauchy dist
    % interesting, for higher dof, this interpolates between Gaussian and
    % Cauchy... for dof=10, it's almost log-linear over 7 standard
    % deviations, and then in flattens out...
    
    % normalize if needed
    m = max(max(obs_prob));
    if m>1
        obs_prob = obs_prob / m;
    end
    
end