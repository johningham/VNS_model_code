function dudt=VNS_vectorise(~,u,p)

% u is the initial state of the populations and needs to be a vector the 
% same length as the number of regions (22).

% p is the parameter set, to allow testing the model with different
% parameters and for doing bifurcation scans. 

    % additional inputs for the thalamic areas
    thal = zeros(22,1);
    thal(3) = -p.TC2RE * (p.a * u(4) + p.b);
    thal(4) = p.RE2RE * (p.a * u(3) + p.b) -p.RE2TC * (p.a * u(4) + p.b);
     
    % change in population outputs over time
    dudt = (p.h - u + p.w * (1 ./ (1+p.epsilon .^ -(u))) + thal) .* p.tau';

end