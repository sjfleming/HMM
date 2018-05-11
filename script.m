% Test the HMM class on nanopore-like simulated data.
% Stephen Fleming
% 4/5/18

clear all
close all

%% set up the sequence

% specify a sequence
%seq = 'RRRRRTTTTTGGGAAATTTTTGGGAAATTTTT';
%seq = 'RRRRRTTTTTGGGAAATTTTTGGGAAATTTTTAGCTGATTACGCTTAGCGGCAATCGGATCCCGATCGATTGGGCTAGCTTAGCGATCGACTGTCTAGGCTTTAAGCGCCCCATTCGAGGTTTAGGATCGGA';
seq = 'TTTTTGGGAAATTTTTGGGAAATTTTTAGCTGATTACGCTTAGCGTTTTCGCAATCGGATCCCGATCGA';

% get the current levels for the sequence
measured_levels = linspace(100,300,numel(seq));
[pA, pA_std, ~, ~, ~] = get_model_levels_M2(seq, measured_levels);

% scale and offset the data from the underlying model
scale = 0.85;
offset = 20;

%% generate simulated data

N = 2000;
t = exprnd(1,1,numel(pA));
t(1) = 2; % necessary minimum condition for indexing to work
t = cumsum(t);
t = round(t/t(end)*(numel(pA)-1)*N/numel(pA));
t = t(1:end-1); % indices where level transitions happen
t = [1, t, N];

d = zeros(N,1);
states = cell(1,numel(pA));
for i = 1:numel(t)-1
    % model states cell array
    states{i}.level_mean = pA(i);
    states{i}.level_stdv = pA_std(i);
    % data
    if t(i+1)<=t(i)
        continue;
    end
    d(t(i):t(i+1)) = pA(i)*scale + offset + randn(t(i+1)-t(i)+1,1)*pA_std(i);
end

% put in noise! (maybe...)
nnoise = 100;
p_noise = 1/N;
if binornd(1,0)
    snoise = max(1,floor((N-nnoise)*rand));
    d(snoise:(snoise+nnoise-1)) = 250*rand + randn(nnoise,1)*2;
end

%% set up the HMM

p_stay = (N-numel(pA))/N;
p_forward = numel(pA)/N;
T = transition_matrix(numel(pA),p_forward^30,p_stay,p_forward,p_forward^10);
fac = p_forward; % down by this factor since you'll get on avg this many emissions per level
hmm = HMM('data',d,'transition',T,'emission',@(x) emission_prob(x,states,p_noise) * fac);

%% do a viterbi alignment

sizeoffont = 14;

%hmm.viterbi;
hmm.EM_viterbi;
figure
plot(d)
ylabel('Current (pA)')
xlabel('Time index')
title('Simulated data')
set(gca,'fontsize',sizeoffont,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0])

figure
plot(hmm.viterbi_alignment.states,'.-')
hold on
stairs(t,[1:numel(pA), numel(pA)],'r')
ylabel('Model state')
xlabel('Time index')
set(gca,'fontsize',sizeoffont,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0])

hmm.plot_alignment('Simulated data after Viterbi alignment');
ylabel('Current (pA)')

display(['Found data scale = ' num2str(1/hmm.data_scaling(1),3) ', and offset = ' ...
    num2str(-hmm.data_scaling(2)/hmm.data_scaling(1),3)])
