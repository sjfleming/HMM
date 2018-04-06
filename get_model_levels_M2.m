function [pA, pA_std, std_mean, scale, offset] = get_model_levels_M2(seq, measured_levels)
% get_model_levels_M2(seq, measured_levels) returns current levels for
% sequence 'seq' that are predicted by M2-MspA measurements by Laszlo et al
% and scaled to mesh with the empirical CDF of the levels passed in,
% 'measured_levels'.
% Stephen Fleming, 6/1/16

    % Takes in a sequence in the form of 'ATCGTTCAAAGC...' with 'R'
    % denoting abasic sites

    % generate "states"
    k = 5;
    states = get_states(nt2int(seq), k); % k for k-mer

    % get the current levels from my own M2 model
    d = load('M2_model.mat');
    pA = nan(size(states));
    std_mean = ones(size(states));
    
    % real parts are the parts without any abasics
    pA(imag(states)==0) = d.M2_model.level_mean(states(imag(states)==0)); % convert into something like pA
    pA_std(imag(states)==0) = d.M2_model.level_stdv(states(imag(states)==0));
    std_mean(imag(states)==0) = d.M2_model.sd_mean(states(imag(states)==0));
    
    % imaginary parts are the abasic parts (we have to search the model)
    [~,iminds] = find(imag(states)~=0);
    for i = iminds
        strstate = seq(i:(i+k-1)); % get the string sequence for abasic kmer
        abasicind = find(cellfun(@(x) strcmp(strstate,x), cellstr(d.M2_model_abasics.kmer')),1,'first');
        pA(i) = d.M2_model_abasics.level_mean(abasicind); % find the matching string and get level info
        pA_std(i) = d.M2_model_abasics.level_stdv(abasicind);
    end
    
%     % scale the current levels to match the scaling of the measured levels
%     
%     [cdf1,lev1] = ecdf(pA);
%     [cdf2,lev2] = ecdf(measured_levels);
%     
%     function value = fun(m,b,cdf1,cdf2,lev1,lev2)
%         value = sum((interp1(cdf2,lev2,linspace(0.05,0.95,10),'linear','extrap')'-(m.*interp1(cdf1,lev1,linspace(0.05,0.95,10),'linear','extrap')'+b)).^2);
%     end
%     
%     params = fminsearch(@(p) fun(p(1),p(2),cdf1,cdf2,lev1,lev2), [std(lev1)/std(lev2), mean(lev1)-mean(lev2)]);
%     scale = params(1);
%     offset = params(2);
%     
%     pA = pA * scale + offset;
%     pA_std = pA_std * scale;
    
    % scale the current levels to match the scaling of the measured levels
    
    xx = 0:400;
    [y1,x] = hist(pA,0:5:400);
    y1 = y1/sum(y1);
    y1 = min(y1,0.03);
    [y2,~] = hist(measured_levels,x);
    y2 = y2/sum(y2);
    y2 = min(y2,0.03);
    y1 = interp1(x,y1,0:400);
    y2 = interp1(x,y2,0:400);
    
    function value = fun(m,b,y1,y2)
        yy1 = interp1((0:400)*m+b,y1,0:400,'linear','extrap'); % scale the levels and reinterpolate
        value = sum((y2-yy1).^2);
    end
    
    params = fminsearch(@(p) fun(p(1),p(2),y1,y2), [range(measured_levels)/range(pA), min(measured_levels)-min(pA)*range(measured_levels)/range(pA)]);
    scale = params(1);
    offset = params(2);
    
    pA = pA * scale + offset;
    pA_std = pA_std * scale;

end