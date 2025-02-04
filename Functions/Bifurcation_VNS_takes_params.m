function fig = Bifurcation_VNS_takes_params(p,param_to_change,paramrange,repeats)
% Bifurcation diagram creator for the deterministic version of our VNS
% model.

% If param_to_change is a scaler integer from 1 to 22, the corresponding 
% population background activity level (h) will automatically be swept. 
% If passed a vector of two such values, will sweepconnection weight (w). 

ignore4time = 6; % sets time after which things get measured

if nargin < 4
    repeats = 1;
end

set(0,'defaultfigurecolor',[1 1 1]) % set figure background to white.

%initial condition near the fixed point
%                      S1_PY, SI_IN,     TC,    RE,INS_PY,INS_IN,ACC_PY,ACC_IN,
near_FP_conditions = [0.1724,0.1787,-0.0818,0.2775,0.0724,0.0787,0.0724,0.0787,...
    0.0724,0.0787,0.0724,0.0787,0.0724,0.0787,0.0724,0.0787,0.0724,0.0787,0.0724,0.0787,0.1724,0.1787];
%   PFC_PY,PFC_IN,Amy_PY,Amy_In,Hyp_Ex,Hyp_In, LC_PY, LC_In,DRN_PY,DRN_In, PB_PY, PB_In,STN_PY,STN_In

init_cond = near_FP_conditions; 

figure; hold on; % set up blank canvas and keep plotting on it.

for rr = 1:repeats % 1 by default, but can loop for multiple start points by including repeats > 1

    % for each forward/back pair remember starting point so the same each way
    reset_cond = init_cond; 
    
    % initialise variables to plot, these can vary in length so starting empty.
    plot_params_fmin = [];
    plot_params_fmax = [];
    plot_params_bmin = [];
    plot_params_bmax = [];
    s1_min_f = [];
    s1_max_f = [];
    s1_min_b = [];
    s1_max_b = [];
        
    %% forward scan:

    % run model once after randomising to stabilise system
    [~,u]=ode45(@(t,u)VNS_vectorise(t,u,p),[0 200],init_cond);%background state
    init_cond = u(end,:);

    for i=1:length(paramrange)
        
        if isscalar(param_to_change) %  if it is a scalar then it is a static input
            p.h(param_to_change) = paramrange(i); % need to work out how to change this in a sensible way
        else % it should be a weight in the form [x,y]
            p.w(param_to_change) = paramrange(i);
        end

        % run the model:
        [t,u]=ode45(@(t,u)VNS_vectorise(t,u,p),[0 8],init_cond);%background state
        
        % save the min/max for each u:
        umin_f(i,:) = min(u(t>ignore4time,:));
        umax_f(i,:) = max(u(t>ignore4time,:));
        % save for s1 (combining PY and IN as in Peter's paper sims)
        s1 = mean([u(:,1),u(:,2)],2);
        s1_min = findpeaks(-s1(t>ignore4time));
        s1_min = -unique(round(s1_min,2));
        % find maximums
        s1_max = findpeaks(s1(t>ignore4time));
        s1_max = unique(round(s1_max,2));
        
        if length(s1_min)~=length(s1_max) % this means something has gone a little squiffy.
            warning(['Uneven min/maxes at parameter =',num2str(paramrange(i)),... 
                ' in forwards run, just taking the pure max/min rather than peak values for safety.'])
            s1_min = min(s1(t>ignore4time));
            s1_max = max(s1(t>ignore4time));
        elseif isempty(s1_min) || isempty(s1_max) % if there are no peaks at all
            s1_min = min(s1(t>ignore4time));
            s1_max = max(s1(t>ignore4time));
        end
                
        s1_min_f = [s1_min_f;s1_min];
        s1_max_f = [s1_max_f;s1_max];
        
        if length(s1_min)>1 % if we have a SWD this will be 2 or more
            % duplicate the current parameter value
            plot_params_fmin = [plot_params_fmin;paramrange(i)*ones(length(s1_min),1)];
        else
            plot_params_fmin = [plot_params_fmin;paramrange(i)];
        end

        if length(s1_max)>1 % if we have a SWD this will be 2 or more
            % duplicate the current parameter value
            plot_params_fmax = [plot_params_fmax;paramrange(i)*ones(length(s1_max),1)];
        else
            plot_params_fmax = [plot_params_fmax;paramrange(i)];
        end
          
        % set new initial conditions based on the output from the last run:
        init_cond = u(end,:);
    end
    
    %% backward scan:
    
    init_cond = reset_cond;
    params_b = fliplr(paramrange); % reverse the parameters

    % run model once after randomising to stabilise system
    [~,u]=ode45(@(t,u)VNS_vectorise(t,u,p),[0 200],init_cond);%background state
    init_cond = u(end,:);

    for i=1:length(params_b)
        
        if isscalar(param_to_change) %  if it is a scalar then it is a static input
            p.h(param_to_change) = params_b(i); % need to work out how to change this in a sensible way
        else % it should be a weight in the form [x,y]
            p.w(param_to_change) = params_b(i);
        end
        
        % run the model:
        [t,u]=ode45(@(t,u)VNS_vectorise(t,u,p),[0 8],init_cond);%background state
        
        % save the min/max for each u:
        umin_b(i,:) = min(u(t>ignore4time,:));
        umax_b(i,:) = max(u(t>ignore4time,:));
        % save for s1 (combining PY and IN as in Peter's paper sims)
        s1 = mean([u(:,1),u(:,2)],2);
        s1_min = findpeaks(-s1(t>ignore4time));
        s1_min = -unique(round(s1_min,2));
        % find maximums
        s1_max = findpeaks(s1(t>ignore4time));
        s1_max = unique(round(s1_max,2));
        
        if length(s1_min)~=length(s1_max) % this means something has gone a little squiffy.
            warning(['Uneven min/maxes at parameter =',num2str(paramrange(i)),...
                ' in backwards run, just taking the pure max/min rather than peak values for safety.'])
            s1_min = min(s1(t>ignore4time));
            s1_max = max(s1(t>ignore4time));
        elseif isempty(s1_min) || isempty(s1_max) % if there are no peaks at all
            s1_min = min(s1(t>ignore4time));
            s1_max = max(s1(t>ignore4time));
        end
        
        s1_min_b = [s1_min_b;s1_min];
        s1_max_b = [s1_max_b;s1_max];

        if length(s1_min)>1 % if we have a SWD this will be 2 or more
            % duplicate the current parameter value
            plot_params_bmin = [plot_params_bmin;params_b(i)*ones(length(s1_min),1)];
        else
            plot_params_bmin = [plot_params_bmin;params_b(i)];
        end

        if length(s1_max)>1 % if we have a SWD this will be 2 or more
            % duplicate the current parameter value
            plot_params_bmax = [plot_params_bmax;params_b(i)*ones(length(s1_max),1)];
        else
            plot_params_bmax = [plot_params_bmax;params_b(i)];
        end
        
        % set new initial conditions based on the output from the last run:
        init_cond = u(end,:);
    end
    
    %% Plotting bifurcation diagram: S1
    hold on
    scatter(plot_params_fmin,s1_min_f,'or','filled','MarkerFaceAlpha',0.1);
    scatter(plot_params_fmax,s1_max_f,'or','filled','MarkerFaceAlpha',0.1);
    scatter(plot_params_bmin,s1_min_b,'ob','filled','MarkerFaceAlpha',0.05);
    scatter(plot_params_bmax,s1_max_b,'ob','filled','MarkerFaceAlpha',0.05);
    if length(param_to_change)>1
        xlabel(['Connection weight ',num2str(param_to_change)])
    else
        xlabel(['Static input ',num2str(param_to_change)]);
    end
    ylabel('Termination Points')
    title('Bifurcation diagram, S1')
    
    % First sweeps run from FP, second sweeps from LC.
    % Subsequent sweeps from random points at either end.
    % (Default value 2, used for paper figure)
    if rr == 1
        init_cond = zeros(1,22);
    else
        init_cond = randn(1,22);
    end
end

fig = gca;

end