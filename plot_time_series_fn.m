function fig = plot_time_series_fn(p, init_cond)
% Takes a set of parameters and a set of initial conditions and returns 
% figure showing a section of time series after equilibrium reached.
    
    % run model once after randomising to stabilise system
    [~,u]=ode45(@(t,u)VNS_vectorise(t,u,p),[0 6],init_cond); 
    init_cond = u(end,:);
    
    % Run the model:
    [t,u]=ode45(@(t,u)VNS_vectorise(t,u,p),[0 2],init_cond);
    
    % Plotting model results - S1 EEG
    % Combine S1 inhibitory and excitatory populations to produce simulated EEG
    s1 = mean([u(:,1),u(:,2)],2);
    figure
    plot(t,s1,'k')
    xlabel('Time (sec)','FontSize',20)
    ylabel('S1','FontSize',20)
    title('VNS vectorised')

    fig = gca;
end