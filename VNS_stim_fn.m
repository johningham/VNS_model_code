function dudt=VNS_stim_fn(t,u,p)

% Copied and adapted from VNSfn_stoch_vec_Euler_stim
% Main change is that p.stimVal is used with existing p.h(21) and (22) 
% values and the current (overall) timestep, frequency, pulse width and
% duty cycle info to generate complex pattern of VNS in keeping with
% the way it is used clinically.

% calculate the timestep over all epochs
masterStep = (p.epoch - 1) * (p.endtime/p.dt) + (t/p.dt);

% Is the stim on relative to pulse frequency?
smallPeriodSteps = round(1/(p.dt * p.freq)); 
stepInSmPeriod = mod(masterStep,smallPeriodSteps);
pwSteps = round(p.pw/(p.dt*10^6));
smallWave = stepInSmPeriod < pwSteps;

% do similar with cycleOn and Off. 
bigPeriodSteps = round((p.cycleOff+p.cycleOn)/p.dt);
stepInBigPeriod = mod(masterStep,bigPeriodSteps);
cyOnSteps = round(p.cycleOn/p.dt);
bigWave = stepInBigPeriod < cyOnSteps;

% Have bigWave switch smallWave
stimOn = bigWave * smallWave;

% Add any stim required to NTS offsets.
p.h(21) = p.h(21) + p.stimVal(1)*stimOn; % Py
p.h(22) = p.h(22) + p.stimVal(2)*stimOn; % Inh

% additional inputs for the thalamic areas
thal = zeros(22,1);
thal(3) = -0.6*(p.s*u(4)+.5);
thal(4) = 10.5*(p.s*u(3)+.5) -0.2*(p.s*u(4)+.5);

% stochastic version (using indexing of passed noise vector):
ix = int32((t+p.dt)/p.dt);
dudt = (p.h - u + p.w * (1./(1+250000.^-(u))) + thal + p.noisevecs(ix,:)').*p.tau';

end
