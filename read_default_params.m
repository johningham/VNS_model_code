function [p] = read_default_params()

    % Sets the default connectivity parameters between regions. 
    % (Weights needs to be a matrix of region x region with entries only 
    %  where there are connections present.)
    struct = load('VNSconnectivity.mat');
    w = struct.mat;

    % note that weights between TC and RE regions are not included above,
    % due to a different set of equations being used. 
    % These are handled separately below...
    p.TC2RE = 0.6;
    p.RE2TC = 0.2;
    p.RE2RE = 10.5;


    % Default static input levels
    h(1) = -.35;      % S1 PY
    h(2) = -3.4;      % S1 IN
    h(3) = -2;        % TC
    h(4) = -12;       % RE 
    h(5:2:21) = h(1); % set all other PYs to the same base input by default.
    h(6:2:22) = h(2); % set all other INs to the same base input by default.
    
    h=h'; % Flip to column vector
    
    a=2.8;       % thalamic linear activation function gradient
    b=0.5;       % thalamic linear activation function y-intercept
    epsilon = 250000; 
    % AffIn = zeros(1,22);
    % AffIn(21) = 1;
    % for the deterministic to replicate VNS Amari this will be all zeros
    % except for NTS excitatory which is set to 1.
        
   % Time scale parameters
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

   
   p.w = w;
   p.h = h;
   p.a = a;
   p.b = b;
   p.epsilon = epsilon;
   % p.AffIn = AffIn;
   p.tau = tau;
