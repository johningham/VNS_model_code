function [p] = read_default_params()
% Allows setup of standard setup of parameters used as a starting point by
% all versions of the model.

    % Sets the default connection weights between populations. 
    % Weights are stored as a 22x22 matrix, w, where w(i,j) represents the
    % weight of a connection from poluation j to population i.
    struct = load('VNSconnectivity.mat');
    w = struct.mat;

    % Note that weights between TC and RE regions are not included above,
    % due to a different set of equations being used. 
    % These are handled separately, with the following values...
    p.TC2RE = 0.6;
    p.RE2TC = 0.2;
    p.RE2RE = 10.5;

    % Default static input levels
    h(1)  = -0.35;    % S1_PY
    h(2)  = -3.4;     % S1_IN
    h(3)  = -2;     % TC
    h(4)  = -12;      % RE 
    h(5:2:21) = h(1); % Set all other PYs to the same base input by default.
    h(6:2:22) = h(2); % Set all other INs to the same base input by default.
        
    % Time scale parameters for each population
    tau(1) = 26;       
    tau(2) = 32.5;    
    tau(3) = 2.6;    
    tau(4) = 2.6;    
    tau(5) = 26;       
    tau(6) = 32.5;    
    tau(7) = 26;       
    tau(8) = 32.5;    
    tau(9) = 26;       
    tau(10) = 32.5;    
    tau(11) = 26;       
    tau(12) = 32.5;    
    tau(13) = 26;       
    tau(14) = 32.5;    
    tau(15) = 26;       
    tau(16) = 32.5;    
    tau(17) = 26;       
    tau(18) = 32.5;    
    tau(19) = 26;       
    tau(20) = 32.5;    
    tau(21) = 26;       
    tau(22) = 32.5;    
    
    a=2.8;               % Thalamic linear activation function gradient.
    b=0.5;               % Thalamic linear activation function y-intercept.
    epsilon = 250000;    % Determines sigmoid activation function gradient.  

    p.w = w;
    p.h = h'; % (flips to column vector)
    p.a = a;
    p.b = b;
    p.epsilon = epsilon;
    p.tau = tau;
    
end