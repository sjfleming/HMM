% Test the HMM class on nanopore-like simulated data.
% Stephen Fleming
% 4/5/18

clear all
close all

%% set up the sequence

% specify a sequence
seq = 'RRRRRTTTTTGGGAAATTTTTGGGAAATTTTT';
%seq = 'RRRRRTTTTTGGGAAATTTTTGGGAAATTTTTAGCTGATTACGCTTAGCGGCAATCGGATCCCGATCGATTGGGCTAGCTTAGCGATCGACTGTCTAGGCTTTAAGCGCCCCATTCGAGGTTTAGGATCGGA';

% get the current levels for the sequence
measured_levels = linspace(100,300,numel(seq));
[pA, pA_std, ~, ~, ~] = get_model_levels_M2(seq, measured_levels);

% scale and offset data from model
scale = 0.85;
offset = 20;

%% generate simulated data

N = 5000;
t = exprnd(1,1,numel(pA));
t = cumsum(t);
t = round(t/t(end)*(numel(pA)-1)*N/numel(pA));
t = t(1:end-1); % indices where level transitions happen
t = [1, t, N];

d = zeros(N,1);
states = cell(1,numel(pA));
for i = 1:numel(t)-1
    % data
    d(t(i):t(i+1)) = pA(i)*scale + offset + randn(t(i+1)-t(i)+1,1)*pA_std(i);
    % model states cell array
    states{i}.level_mean = pA(i);
    states{i}.level_stdv = pA_std(i);
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
T = transition_matrix(numel(pA),p_forward^5,p_stay,p_forward,p_forward^2);
fac = p_forward; % down by this factor since you'll get on avg this many emissions per level
hmm = HMM('data',d,'transition',T,'emission',@(x) emission_prob(x,states,p_noise) * fac);

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
plot(hmm.viterbi_alignment.states,'.-')
hold on
stairs(t,[1:numel(pA), numel(pA)],'r')
ylabel('Model state')
xlabel('Time index')
set(gca,'fontsize',sizeoffont,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0])

figure
a = 0;
threshold = log10(p_noise)+0.01;
noise_logic = hmm.viterbi_alignment.log_emissions < threshold;
xx = 1:numel(d);
c = colormap(lines(numel(pA)));
for i = 1:numel(pA)
    state_logic = hmm.viterbi_alignment.states==i;
    plot(xx(state_logic & ~noise_logic),d(state_logic & ~noise_logic),'color',c(i,:))
    hold on
    plot(xx(state_logic & noise_logic),d(state_logic & noise_logic),'o','color',c(i,:))
    if sum(hmm.viterbi_alignment.states==i)>0
        text(a+sum(hmm.viterbi_alignment.states==i)/5,states{i}.level_mean+15,num2str(i),'fontsize',9)
    end
    a = a+sum(hmm.viterbi_alignment.states==i);
end
ylabel('Current (pA)')
xlabel('Time index')
title('Simulated data after Viterbi alignment')
set(gca,'fontsize',sizeoffont,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0])

