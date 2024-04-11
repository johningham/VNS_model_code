function [p] = read_default_params()


    % default connectivity parameters
    
    % weights needs to be a matrix of regionxregion with entries only where
    % there are connections present
    struct = load('VNSconnectivity.mat');
    w = struct.mat;

    % default static input levels
    h(1) = -.35;      % S1 PY
    h(2) = -3.4;      % S1 IN
    h(3) = -2.05;     % TC, -2 for bistable for ode45; -2.2 for excitable for ode45. -2.05 for full model
    h(4) = -12;       % RE -5 in default model, -12 for full model
    h(5:2:21) = h(1); % set all PYs to the same base input by default for now.
    h(6:2:22) = h(2); % set all INs to the same base input by default for now.
    
    h=h'; % flip to column vector
    
    s=2.8;     % thalamic input multiplier
    a=1;       % thalamic tau (time delay) multiplier
    AffIn = zeros(1,22);
    AffIn(21) = 1;
    % repmat([0.1,0],1,11); % Afferent input scaler,
    % for the noisy sim setting the default to 0.1* input to excitatory populations and 0*input to inhibitory
    % for the deterministic to replicate VNS Amari this will be all zeros
    % except for NTS excitatory which is set to 1.
        
   % Time scale parameters
   tau(1)=1*26;       %
   tau(2)=1.25*26;    %
   tau(3)=.1*a*26;    %
   tau(4)=.1*a*26;    %
   tau(5)=1*26;       %
   tau(6)=1.25*26;    %
   tau(7)=1*26;       %
   tau(8)=1.25*26;    %
   tau(9)=1*26;       %
   tau(10)=1.25*26;    %
   tau(11)=1*26;       %
   tau(12)=1.25*26;    %
   tau(13)=1*26;       %
   tau(14)=1.25*26;    %
   tau(15)=1*26;       %
   tau(16)=1.25*26;    %
   tau(17)=1*26;       %
   tau(18)=1.25*26;    %
   tau(19)=1*26;       %
   tau(20)=1.25*26;    %
   tau(21)=1*26;       %
   tau(22)=1.25*26;    %
   
   p.w = w;
   p.h = h;
   p.s = s;
   p.AffIn = AffIn;
   p.tau = tau;
