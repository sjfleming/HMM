classdef HMM < handle
% HMM
% Class that handles common computations for hidden Markov models.
% Required inputs: (name, value pairs)
% 'transition': a square matrix of transition probabilities between states
% in the model.  matrix is N by N.
% 'emission': a function handle that, when evaluated on a single data
% values, returns an N-element array of emission probabilities of that
% observation from each state of the model.
% 'data': a time-series array of M data points.
    
    properties
        % these properties go with the HMM object
        
        % data and model
        % sigdata = []; % SignalData object associated with data file, or raw data
        data = [];
        T = []; % transition matrix
        E = @(x) []; % emission probability generating function
        init = []; % initial state occupancy probabilities
        data_scaling = [1,0];
        
        % generated quantities
        viterbi_alignment = [];
        
    end
    
    methods % can be called directly by the user on a particular analysis object
        
        % constructor
        function obj = HMM(varargin)
            % handle inputs
            in = obj.parse_inputs(varargin{:});
            
            % initialize the HMM object by filling in property values
            obj.data = in.data;
            obj.T = in.transition;
            obj.E = in.emission;
            obj.init = in.starting_prob;
            
        end
        
        function alignment = viterbi(obj)
            % use the Viterbi algorithm to calculate the best possible
            % alignment of the data to the model states, given the
            % transition matrix
            
            % use default initial probs if none given
            if isempty(obj.init)
                n = size(obj.T,1);
                initial = exp(-(1:n)/5e-1)'/sum(exp(-(1:n)/5e-1));
                if numel(initial)>n
                    initial = initial(1:n) / sum(initial(1:n)); % concatenate and fix to sum to 1
                end
            else
                initial = obj.init;
            end
            loginit = log10(initial);
            
            % initialize the big matrix 'dp' to trace path
            % rows are all possible states in the model
            % columns are the observations
            % dimension 3 contains: prob, pointer i
            dp = zeros(size(obj.T,1), numel(obj.data), 2);
            % fill in first column, for first observation
            dp(:,1,1) = loginit + log10(obj.E(obj.data(1)));

            % fill in the big matrix
            % each element (i,j,1) stores the log prob of the most likely path so far
            % and a pointer (i,j,2) to the previous row index of that path
            logT = log10(obj.T);
            for j = 2:numel(obj.data) % columns are observations

                % calculate the emission probs for every model state
                % given this observation (a column vector)
                logEm = log10(obj.E(obj.data(j)));
                
                for i = 1:size(obj.T,1) % rows are all possible states

                    % all possible previous states times transition to this state
                    % times emission for this state
                    [m, ind] = max( dp(:,j-1,1) + logT(:,i) + logEm(i)); % sum log probabilities
                    %m = m + logEm(i); % also take into account the emission prob
                    dp(i,j,1) = m; % the probability of the maximally probable path to here
                    dp(i,j,2) = ind; % row index

                end

            end
            
            % get the most probable path by finding the last most probable state
            [~,ind] = max(dp(:,end,1)); % state index
            z = nan(1,numel(obj.data)); % best path state indices
            z(end) = ind;

            % trace back through the big matrix to get the sequence of states
            for j = numel(obj.data):-1:2

                z(j-1) = dp(z(j),j,2); % pointer to the previous row index, i.e. state index

            end
            
            obj.viterbi_alignment = z;
            alignment = z;
            
        end
        
        function plot_model(obj)
            % plot the states of the HMM as a directed graphical model,
            % with edges weighted by their transition probabilities
            
            % create a directed graph object, using Matlab's built-ins
            % where edges are the transition matrix elements
            Tr = obj.T;
            Tr(Tr<1e-4) = 0; % set small probability transitions to zero
            ind = round(size(Tr,1)/2); % row index near the middle
            inds = 1:size(Tr,1);
            Tr2 = Tr;
            Tr2(inds~=ind,:) = 0; % select only that one middle row
            
            % plot full graph, with edges having thickness proportional to prob
            figure
            g = digraph(Tr);
            f = plot(g);
            f.NodeColor = [0 0 0];
            f.MarkerSize = 8;
            g.Edges.LWidths = 7*g.Edges.Weight/max(g.Edges.Weight);
            f.LineWidth = g.Edges.LWidths;
            title('Full model of state transitions')
            
            % plot full graph, with edges having thickness proportional to prob
            figure
            g2 = digraph(Tr2);
            f2 = plot(g2);
            f2.NodeColor = [0 0 0];
            f2.MarkerSize = 8;
            g2.Edges.LWidths = 7*g2.Edges.Weight/max(g2.Edges.Weight);
            f2.LineWidth = g2.Edges.LWidths;
            title(['Transitions from state ' num2str(ind) ' in model'])
            
        end
        
    end
    
    methods (Access = public) % only called by class methods
        
        
        
        function [scale, offset] = data_to_model_transformation(obj)
            % suppose there's a linear function mapping the "true" model to
            % the model we're using (scale and offset parameters) that are
            % unknown.  the data were taken at a different temperature,
            % salt concentration, whatever, so the scaling is off.  we
            % don't have explicit access to the model state information,
            % only the emission probability function.  we would like to
            % apply a scaling and offset to the model states, but we can
            % just as easily apply the inverse to the data.  then we'll
            % know that in reality, we need to do the inverse scaling to
            % the model.
            % this does an L1 regression to fit the data to the model.
            
            % get the max likelihood emission values for each model state
            model_means = obj.get_model_values();
            
            % if we haven't aligned yet, make an inital guess about the
            % model scaling and offset by fitting a line
            if isempty(obj.viterbi_alignment)
            
                % randomly grab points from model and data
                m = datasample(model_means',500);
                d = datasample(obj.data,500);

                % sort them in order, adding noise to the model
                blur = range(obj.data)/size(obj.T,1)/2; % empirical estimate
                m = sort(m + randn(size(m))*blur);
                d = sort(d);
            
            else
                % we do have a viterbi alignment, so use that information
                % to do a linear fit of model to data
                
                m = sort(model_means)';
                
                
            end
            
            % fit a line and return the best fit params
            %[scale, offset] = polyfit(m,d,1); % polynomial order 1, a line
            p = robustfit(m,d); % better at ignoring outliers
            offset = p(1);
            scale = p(2);
            
        end
        
        function model_means = get_model_values(obj)
            % we don't know the model explicitly...
            % approximate the model levels by doing a max over emissions
            
            current = (min(obj.data)-100):(max(obj.data)+100); % currents
            logEm = zeros(numel(current),size(obj.T,1)); % allocate space
            for j = 1:numel(current) % for each current in the data range
                logEm(j,:) = log10(obj.E(current(j))); % get emissions
            end
            [~,inds] = max(logEm,[],1); % max prob emission index
            model_means = current(inds); % max prob emission current
            
        end
        
        function in = parse_inputs(obj, varargin)
            % parse all inputs so all methods can use them easily
            p = inputParser;
            
            % defaults and checks
            check_trans_matrix = @(x) all([isnumeric(x), ...
                all(x>=0), ...
                all(x<=1), ...
                size(x,1)==size(x,2)]);
            check_emission = @(x) all([isnumeric(feval(x,[1,2,3])), ...
                all(feval(x,[1,2,3])>=0), ...
                all(feval(x,[1,2,3])<=1)]);%, ...
                %numel(feval(x,[1,2,3]))==3]);
            check_init = @(x) all([isnumeric(x), ...
                all(x>=0)]);
            
            % set up the inputs
            addOptional(p, 'data', [], @isnumeric); % data observed
            addOptional(p, 'emission', [], check_emission); % emission prob function handle
            addOptional(p, 'transition', [], check_trans_matrix); % transition prob matrix
            addOptional(p, 'starting_prob', [], check_init); % starting probability array
            %addOptional(p, 'states',
            
            % parse
            parse(p,varargin{:});
            in = p.Results;
            
            % impose normalization on transition probabilities
            % i.e. rows must sum to one
            in.transition = obj.row_normalize(in.transition);
            in.starting_prob = in.starting_prob / sum(in.starting_prob); % normalize
            
        end
        
        function n = row_normalize(~, m)
            % normalize rows of matrix m, so they add to one
            % return normalized matrix
            n = m ./ repmat(sum(m,2), 1, size(m,1));
        end
        
    end
    
end