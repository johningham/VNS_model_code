function [t,u] = vectorised_eulersolver(funcname,u_init,dt,tend)

% Simple solver using Euler-Maruyama for our vectorised code to run the
% stochastic versions of the model.
    
    % The time vector is in steps of dt up to tend 
    t = 0:dt:tend;
    
    % initialise the solution u, parametersxtime
    u = zeros(length(u_init),length(t));
    u(:,1) = u_init;
    
    % euler solver
    for i = 1:length(t)-1
        
        dudt = funcname(t(i),u(:,i));
        u(:,i+1) = u(:,i) + dudt*dt;
        
    end   

    u = u'; % transpose to match the ode output format for u.

end

