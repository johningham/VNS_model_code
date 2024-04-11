function dudt=VNS_vectorise(t,u,p)
% For this version of the VNS function, u should be the initial conditions 
% which needs to be a vector the same length as the number of regions (22).
% p is the parameter to set, to allow testing the model with different
% parameters and for doing bifurcation scans. 
% If p is empty the default parameters will be used.
% Otherwise, p should be a structure with the fields:
% p.h should be a boolean, 0 means no inputs need changing, 1 means change a value.
% p.w should be a boolean, 0 means no inputs need changing, 1 means change a value. 
% p.AffIn should be a boolean, 0 means no inputs need changing, 1 means change a value.
% p.num_h is the input parameter number to set. This can be
% more than one, but if so there must be a corresponding number of p.val_h.
% p.num_w is the corresponding weight parameter number to set. This can be
% more than one, but if so there must be a corresponding number of p.val_h.
% p.val_h is the value to set for the changed input parameter. Should be
% either a scalar or a vector with the same number of inputs as p.num_h.
% p.val_w is the value to set for the changed weight parameter. Should be
% either a scalar or a vector with the same number of inputs as p.num_w.
% p.val_AffIn is the value to set for the changed afferent input. Should be
% a scalar.
% There is scope to extend this to vary taus, s or other parameters.
%
% This version assumes that p will be passed in, and no internal defaults
% are provided. For default parameters the readDefaultVNSParams function
% can be called and used to create the p input.
%
% This is an adaptation of VNS_vectorise to keep the clean vectorised
% computation but adding a noise term. The default AffIn from
% readDefaultVNSParams will add this noise only to the excitatory
% populations.
%
% FT 2023


% additional inputs for the thalamic areas
thal = zeros(22,1);
thal(3) = -0.6*(p.s*u(4)+.5);
thal(4) = 10.5*(p.s*u(3)+.5) -0.2*(p.s*u(4)+.5);

% change in population outputs over time
dudt = (p.h - u + p.w * (1./(1+250000.^-(u))) + thal + p.AffIn').*p.tau';

end