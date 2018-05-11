classdef HMM < handle
% HMM
%
% Class that handles common computations for hidden Markov models.
% Required inputs: (name, value pairs)
% 'transition': a square matrix of transition probabilities between states
% in the model.  matrix is N by N.
% 'emission': a function handle that, when evaluated on a single data
% values, returns an N-element array of emission probabilities of that
% observation from each state of the model.
% 'data': a time-series array of M data points.
%
% Built for data sets with tons of observations compared to the number of
% hidden states.  Good for time-series data where there are many noisy
% observations of each hidden state, for example.
%
% Stephen Fleming
% 5/9/2018
    
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
            
            % precompute specified emission probs and store in matrix:
            % this prevents tons of function calls and serves as a look up
            % table, which is necessary for large datasets
            % (limit the table to about 1 million points, for memory)
            num = floor(1e6/size(obj.T,1)); % measurement resolution for emission calculation
            range = [min(obj.scaled_data), max(obj.scaled_data)];
            currents = linspace(range(1),range(2),num);
            div = currents(2) - currents(1);
            logEm = log10(obj.E(currents));

            % fill in the big matrix
            % each element (i,j,1) stores the log prob of the most likely path so far
            % and a pointer (i,j,2) to the previous row index of that path
            logT = log10(obj.T);
            for j = 2:numel(obj.data) % columns are observations

                % calculate the emission probs for (nearly) this current
                data_ind = floor((obj.scaled_data(j) - range(1))/div) + 1;
                
                for i = 1:size(obj.T,1) % rows are all possible states

                    % all possible previous states times transition to this state
                    % times emission for this state
                    [m, ind] = max( dp(:,j-1,1) + logT(:,i) + logEm(i,data_ind) ); % sum log probabilities
                    dp(i,j,1) = m; % the probability of the maximally probable path to here
                    dp(i,j,2) = ind; % row index

                end

            end
            
            % get the most probable path by finding the last most probable state
            [~,ind] = max(dp(:,end,1)); % state index
            z = nan(numel(obj.data),1); % best path state indices
            log_emissions = nan(numel(obj.data),1); % best path emission probs
            z(end) = ind;
            data_ind = floor((obj.scaled_data(j) - range(1))/div) + 1;
            log_emissions(end) = logEm(z(end),data_ind);

            % trace back through the big matrix to get the sequence of states
            for j = numel(obj.data):-1:2

                z(j-1) = dp(z(j),j,2); % pointer to the previous row index, i.e. state index
                
                % trace through the emissions and record the emission probs
                data_ind = floor((obj.scaled_data(j-1) - range(1))/div) + 1;
                log_emissions(j-1) = logEm(z(j-1),data_ind);

            end
            
            obj.viterbi_alignment.states = z;
            obj.viterbi_alignment.log_emissions = log_emissions;
            alignment = obj.viterbi_alignment;
            
        end
        
        function alignment = EM_viterbi(obj, max_iter, tol)
            % a scaling and offset between the model and data is unknown.
            % iteratively (1) estimate scaling and offset, and (2)
            % use the Viterbi algorithm to calculate the best possible
            % alignment of the scaled data to the model states
            
            if nargin < 2
                max_iter = 20;
            end
            if nargin < 3
                tol = 0.0001;
            end
            
            f = figure;
            obj.data_scaling = obj.data_to_model_transformation(true);
            
            for i = 1:max_iter
                display(['EM Viterbi alignment: iteration ' num2str(i)])
                
                obj.viterbi;
                f = obj.plot_alignment(['Iteration ' num2str(i)], f);
                drawnow;
                
                prev = obj.data_scaling;
                obj.data_scaling = obj.data_to_model_transformation;
                delta = max((obj.data_scaling - prev)./prev);
                if delta < tol
                    display('Tolerance achieved')
                    break;
                end
            end
            
            alignment = obj.viterbi_alignment;
            
        end
        
        function f = plot_alignment(obj, plot_title, fig)
            % plot alignment of data to model levels
            
            if isempty(obj.viterbi_alignment)
                display('Alignment has not yet been created.')
                return;
            end
            
            model_means = obj.get_model_values;
            
            if nargin<3
                f = figure;
            else
                f = fig;
                clf(f)
            end
            a = 0;
            threshold = -6; % for "noise" alignment emission probability
            noise_logic = obj.viterbi_alignment.log_emissions < threshold;
            xx = 1:numel(obj.data);
            c = colormap(lines(size(obj.T,1)));
            for i = 1:size(obj.T,1)
                state_logic = obj.viterbi_alignment.states==i;
                plot(xx(state_logic & ~noise_logic),obj.data(state_logic & ~noise_logic),'color',c(i,:))
                hold on
                plot(xx(state_logic & noise_logic),obj.data(state_logic & noise_logic),'o','color',c(i,:))
                if sum(obj.viterbi_alignment.states==i)>0
                    text(a+sum(obj.viterbi_alignment.states==i)/5, ...
                        (model_means(i)-obj.data_scaling(2))/obj.data_scaling(1)+15, num2str(i), 'fontsize',9)
                end
                a = a+sum(obj.viterbi_alignment.states==i);
            end
            ylabel('Data value')
            xlabel('Time index')
            if nargin<2
                title('Viterbi alignment')
            else
                title(plot_title)
            end
            set(gca,'fontsize',14,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0])
            box on;
            
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
        
        function p = data_to_model_transformation(obj, naive)
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
            
            % is this run supposed to be from a naive starting state?
            if nargin<2
                naive = false;
            end
            
            % if we haven't aligned yet, make an inital guess about the
            % model scaling and offset by fitting a line
            if naive || isempty(obj.viterbi_alignment)
            
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
                
                m = model_means(unique(obj.viterbi_alignment.states))';
                d = accumarray(obj.viterbi_alignment.states,obj.data,[],@median);
                d = d(d~=0);
                
            end
            
            % fit a line and return the best fit params
            %[scale, offset] = polyfit(m,d,1); % polynomial order 1, a line
            try
                p = robustfit(d,m); % better at ignoring outliers
                p = fliplr(p'); % [scale, offset]
            catch ex
                p = [m\d, 0];
            end
            %offset = p(2);
            %scale = p(1);
            
        end
        
        function d = scaled_data(obj, inds)
            % return a scaled version of the data according to the object
            % parameters in data_scaling
            
            if nargin < 2
                d = obj.data*obj.data_scaling(1) + obj.data_scaling(2);
            else
                d = obj.data(inds)*obj.data_scaling(1) + obj.data_scaling(2);
            end
            
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