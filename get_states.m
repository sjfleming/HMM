function seqstates = get_states(seq, k)
% Get states from a sequence, using given kmer length (assumed 5).
% 'seq' should be nt2int integers.
% using 'R' for abasics in the nucleotide string works to give the integer
% 5, which is interpreted as abasic.
% If there are abasics, the imaginary part is the number of abasics at
% that kmer location.

%     if (nargin < 2)
%         k = 5;
%     end
    
    % edit 2/8/16, SJF
    % if any of the integers in 'seq' are 5, then this stands for an abasic
    % find these and deal with them first
    abasic_inds = (seq==5);
    if sum(abasic_inds)>0 % there are abasics
        seq_abasics_imaginary = double(seq) - (seq==5)*4;
        seqstates = conv(double(seq_abasics_imaginary)-1,4.^(0:k-1),'valid')+1;
        seqstates = seqstates + conv(double(abasic_inds),ones(1,k)*sqrt(-1),'valid');
        % imaginary number denotes number of abasics at each location
    else
        % usual code
        seqstates = conv(double(seq)-1,4.^(0:k-1),'valid')+1;
    end
end
