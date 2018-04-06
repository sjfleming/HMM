function A = transition_matrix(N, p_back, p_stay, p_forward, p_skip)
% transition_matrix(N, p_back, p_stay, p_skip) returns a transition matrix 
% for an HMM with N states containing the transition probabilities A_ij 
% between state i and state j, computed based on the input probilities of a
% back step, a stay, a forward step, and a skip.
% p_back + p_stay + p_forward + p_skip = 1
% Stephen Fleming, 7/23/16

    % input check
    if p_back + p_stay + p_forward + p_skip ~= 1
        tot = sum([p_back, p_stay, p_forward, p_skip]);
        p_back = p_back / tot;
        p_stay = p_stay / tot;
        p_forward = p_forward / tot;
        p_skip = p_skip / tot;
    end
    
    % transition matrix A, probabilities
    % (probability of going from one state to another)
    block = zeros(N,N);
    for i = 1:N
        row = zeros(1,N);
        % stay
        row(i) = p_stay;
        if i<N
            row(i+1) = p_forward;
        end
        % skip
        skipinds = [i+2, N];
        if skipinds(1)<=N
            skipinds = min(N,skipinds);
            skipinds = skipinds(1):skipinds(2);
            skip = (p_skip.^(1:numel(skipinds)))/sum(p_skip.^(1:numel(skipinds)))*p_skip;
            row(skipinds) = skip;
        end
        % back
        backinds = [1, i-1];
        if backinds(2)>=1
            backinds = max(1,backinds);
            backinds = backinds(1):backinds(2);
            back = (fliplr(p_back.^(1:numel(backinds))))/sum(p_back.^(1:numel(backinds)))*p_back;
            row(backinds) = back;
        end
        
        % ensure A is normalized even though it should be close
        block(i,1:N) = row/(sum(row));
    end
    % A(i,j) is the probability of a transition from state i to state j
    A = block;

end