% Test the HMM class on nanopore-like simulated data.
% Stephen Fleming
% 4/5/18

%% set up the sequence

% specify a sequence
seq = 'RRRRRTTTTTGGGAAATTTTTGGGAAATTTTT';

% get the current levels for the sequence
measured_levels = linspace(100,300,numel(seq));
[pA, pA_std, ~, ~, ~] = get_model_levels_M2(seq, measured_levels);

%% generate simulated data

N = 1000;
t = exprnd(1,1,numel(pA));
t = cumsum(t);
t = round(t/t(end)*(numel(pA)-1)*N/numel(pA));
t = t(1:end-1); % indices where level transitions happen
t = [1, t, N];

d = zeros(N,1);
states = cell(1,numel(pA));
for i = 1:numel(t)-1
    % data
    d(t(i):t(i+1)) = pA(i) + randn(t(i+1)-t(i)+1,1)*pA_std(i);
    % model states cell array
    states{i}.level_mean = pA(i);
    states{i}.level_stdv = pA_std(i);
end

%% set up the HMM

T = transition_matrix(numel(pA),1e-6,(N-numel(pA))/N,1e-2,1e-4);
hmm = HMM('data',d,'transition',T,'emission',@(x) emission_prob(x,states,1e-6));

%% do a viterbi alignment

sizeoffont = 14;

hmm.viterbi;
figure
plot(d)
ylabel('Current (pA)')
xlabel('Time index')
title('Simulated data')
set(gca,'fontsize',sizeoffont,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0])

figure
plot(hmm.viterbi_alignment,'.-')
ylabel('Model state')
xlabel('Time index')
set(gca,'fontsize',sizeoffont,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0])

figure
a = 0;
for i = 1:numel(pA)
    plot(a+1:(a+sum(hmm.viterbi_alignment==i)),d(hmm.viterbi_alignment==i))
    if sum(hmm.viterbi_alignment==i)>0
        text(a+sum(hmm.viterbi_alignment==i)/5,states{i}.level_mean+15,num2str(i),'fontsize',9)
    end
    hold on
    a = a+sum(hmm.viterbi_alignment==i);
end
ylabel('Current (pA)')
xlabel('Time index')
title('Simulated data after Viterbi alignment')
set(gca,'fontsize',sizeoffont,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0])

