function dudt=VNS_vectorise(t,u,p)

% u should be the initial conditions 
% and needs to be a vector the same length as the number of regions (22).

% p is the parameter set, to allow testing the model with different
% parameters and for doing bifurcation scans. 





    % additional inputs for the thalamic areas
    thal = zeros(22,1);
    thal(3) = -p.TC2RE*(p.s*u(4)+.5);
    thal(4) = p.RE2RE*(p.s*u(3)+.5) -p.RE2TC*(p.s*u(4)+.5);
    
    % change in population outputs over time
    dudt = (p.h - u + p.w * (1./(1+250000.^-(u))) + thal).*p.tau';

end