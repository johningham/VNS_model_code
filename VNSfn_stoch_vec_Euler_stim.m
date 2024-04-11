function dudt=VNSfn_stoch_vec_Euler_stim(t,u,p)

% Copied and adapted from VNSfn_stoch_vec_EulerCompatible 4/4/23
% Main change is that this takes p.stimVal which can be used to modify NTS
% input between runs.

% Add any stim required to NTS offsets.
p.h(21) = p.h(21) + p.stimVal(1); % Py
p.h(22) = p.h(22) + p.stimVal(2); % Inh

% additional inputs for the thalamic areas
thal = zeros(22,1);
thal(3) = -0.6*(p.s*u(4)+.5);
thal(4) = 10.5*(p.s*u(3)+.5) -0.2*(p.s*u(4)+.5);

% stochastic version (using indexing of passed noise vector):
ix = int32((t+p.dt)/p.dt);
dudt = (p.h - u + p.w * (1./(1+250000.^-(u))) + thal + p.noisevecs(ix,:)').*p.tau';

end