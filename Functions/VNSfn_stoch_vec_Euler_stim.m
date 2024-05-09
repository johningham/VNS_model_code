function dudt=VNSfn_stoch_vec_Euler_stim(t,u,p)
% p.stimVal can be used to modify NTS input between runs.

    % Add any stim required to NTS offsets.
    p.h(21) = p.h(21) + p.stimVal(1); % Py
    p.h(22) = p.h(22) + p.stimVal(2); % Inh
    
    % additional inputs for the thalamic areas
    thal = zeros(22,1);
    thal(3) = -p.TC2RE * (p.a * u(4) + p.b);
    thal(4) = p.RE2RE * (p.a * u(3) + p.b) - p.RE2TC * (p.a * u(4) + p.b);
    
    % stochastic version (using indexing of passed noise vector):
    ix = int32((t + p.dt) / p.dt);
    dudt = (p.h - u + p.w * (1 ./ (1 + p.epsilon .^ -(u))) + thal ...
        + p.noisevecs(ix,:)') .* p.tau';

end

