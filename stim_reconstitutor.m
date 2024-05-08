% Takes an output file from 'stim_chunker.m' and reconstitutes parts of the
% time series that may be of interest for a specified start point, 
% duration, and combination of parameters. Figure 3 in the paper was also
% made using this script, by setting the VNS stimulation amplitude to zero.
% 
% The code outputs not only time series of the mean S1 populations' time
% series, but also the euclidean distance from the FP, and the smoothed
% version of this used in the seizure detector. This was used in the
% supplementary figure that explained the seizure detection process.
% The code also outputs a figure of time series for the remaining regions
% over the same period with the same parameters.

close all
clear

load('VNS_stim_output_1147915.mat') % data for use in final paper

% settings for plot and data saving
savePlots = false;
saveData = false;

% Choose from available param combinations and specify start time and
% duration:
Oix = 2; 
% (index of parameter One - excitatory component of VNS - low value of 0.2)
Tix = 1;
% (index of parameter Two - inhibitory. Zero is the only value in this set)
Nix = 1;
% (index of the noise scaler value - only one used in this instance: 0.72)
startStep =  245628404;
durStep = 433520;

% Set number of timesteps to reconstruct and plot before and after the
% interval of interest.
startmargin = 50000; % (needs to be no longer than pre-calc timesteps)
endmargin = 100000;

% Fix Y-axis limits
topEdge = 0.5;
bottomEdge = -0.1;

endStep = startStep + durStep;

% unpack p
threshold = p.threshold;
dt = p.dt;
endtime = p.endtime;

stimExVal = p.stimExVals(Oix);
stimInVal = p.stimInVals(Tix);
noiseScaler = p.noiseScalers(Nix);

windsecs = p.windowSecs;
halfwindow = floor(windsecs/(2*dt));
wyndoh = halfwindow*2 +1; % (easier if consistently an odd integer)
clear('nStps')
nStps = cast(p.nSteps,'double');

leftEdge = startStep - startmargin; % (x-limits of plot)
rightEdge = endStep+endmargin;

verystart = cast((startStep - startmargin - halfwindow),'double');
veryend = cast((startStep + durStep + endmargin + halfwindow), 'double');

vsepoch = cast(floor(verystart/nStps),'double') + 1; 
% (epochs numbered from one) 

vslocaltimestep = verystart - ((vsepoch-1) * nStps); 

veepoch = cast((floor(veryend/nStps) + 1), 'double'); 
velocaltimestep = veryend - ((veepoch-1) * nStps); 

initsOTNER = p.init_states; % (paramOne, paramTwo, Noise, Epoch, Region)
initsER = permute((initsOTNER(Oix,Tix,Nix,:,:)),[4 5 1 2 3]);

%% Find precise FP conditions
% (in the absence of noise or stimulation). 

% We have considered doing array of FP conditions for all parameter
% combinations in the sweep. This may have made seizure
% detection more accurate, but had some complications. However,
% the FP remains in a very similar position throughout the range. For the
% purposes of standardisation FP and LC are determined for h values of -1.5
% for NTS_py and -3.4 for NTS_inh, set here:
p.h(21) = -1.5; 

nSteps = int32(10/p.dt + 1);  % (runs for 10 secs to reach equilibrium)
p.stimVal = [0,0]; % (not used, but function needs these passing)
p.noisevecs = zeros(nSteps, 22); 
% (also not used, but function needs these passing too)

nearFP = [0.1724,0.1787,-0.0818,0.2775,0.0724,0.0787,0.0724,0.0787,...
    0.0724,0.0787,0.0724,0.0787,0.0724,0.0787,0.0724,0.0787,0.0724,...
    0.0787,0.0724,0.0787,0.1724,0.1787];
% (a starting state that reliably progresses to the fixed point - across a 
% wide parameter range)

[~,uFP] = vectorised_eulersolver ...
    (@(t,uFP)VNSfn_stoch_vec_Euler_stim(t,uFP,p), nearFP, dt, 10);
exactFP = uFP(end,:);

%% Remaining Setup

% initiate long time series
compositeSeries = zeros(0,22);

% Handle calculations that need to start before time point zero
% (use FP calculating steps - datum keeps track of start point)
if vsepoch == 0
    compositeSeries = cat(1,compositeSeries,uFP); 
    vsepoch = 1;
    datum = -nStps;
else 
    datum = (vsepoch-1) * p.nSteps;
end

% Handle calculations that require exceeding the length of the
% nominal time series. 
% Later will caculate a bit further in all cases 
if veepoch > size(initsER,1)
    veepoch = size(initsER,1);
end

% set back to what was used for main epochs 
nStps = p.nSteps;

% Switch to required constant background NTS levels
p.h(21) = p.baseEx;
p.h(22) = p.baseIn;

% Set the required stimulation level for plot
p.stimVal = [stimExVal stimInVal];

for seed = vsepoch:veepoch
    p.epoch = seed;
    rng(seed, 'twister')
    p.noisevecs = randn(nStps, 22) .* p.noiseCHOICE .* noiseScaler; 
    % (overwrites previous value)
    
    % get the relevant initial conditions to start the epoch
    this_init = initsER((seed),:);

    % send to solver
    [~,u] = vectorised_eulersolver(@(t,u)VNS_stim_fn(t,u,p), this_init, ... 
        dt, endtime-dt);

    % append result to long time series
    compositeSeries = cat(1,compositeSeries,u); 
end

% Calculate a little further on for safety...
seed = seed + 1;
p.epoch = seed;
rng(seed, 'twister')

p.noisevecs = randn(nStps, 22) .* p.noiseCHOICE .* noiseScaler; 
% (overwrites previous value)

% get the relevant initial conditions to start the epoch
this_init = compositeSeries(end,:);

% send to solver
[~,u] = vectorised_eulersolver(@(t,u)VNSfn_stoch_vec_Euler_stim(t,u,p), ...
    this_init, dt, endtime-dt);

% append result to long time series
compositeSeries = cat(1,compositeSeries,u); 

% Calculate LFP (mean of S1 populations)
LFP = mean(compositeSeries(:,1:2),2);

% calculate euclidean distance (using freestanding function)
eucs = eucliser(compositeSeries, p.eucWeights, exactFP);

% get moving averages
smooth_eucs = movmean(eucs,wyndoh);

% chop everything to length
seriesSteps = leftEdge:rightEdge;
preChop = leftEdge - datum;
postChop = (nStps * veepoch) + p.nSteps - rightEdge;
LFP = LFP(preChop:end-postChop);
eucs = eucs(preChop:end-postChop);
smooth_eucs = smooth_eucs(preChop:end-postChop);
compositeSeries = compositeSeries(preChop:end-postChop,:);
startSecs = startStep * dt;
endSecs = endStep * dt;

% plot S1 (included threshold, start of noise if relevant etc)
figure(1)
hold on 
plot(seriesSteps, LFP)
plot(seriesSteps, eucs)
plot(seriesSteps, smooth_eucs,'m')
yline(threshold)
xline(startStep)
xline(endStep)
xlim([leftEdge rightEdge])
xticks([startStep, endStep])
xticklabels([startSecs, endSecs])
xlabel('time (s)')
ylim([bottomEdge topEdge])
f1 = gcf;
f1.Position = [200 200 2000 700];

%% Plotting other regions
stepsToPlot = size(seriesSteps,2);
LFP_Ex = zeros(stepsToPlot, 11);
LFP_In = zeros(stepsToPlot, 11);
count = 1;

for i = 1:2:21
    LFP_Ex(:,count) = compositeSeries(:,i);
    LFP_In(:,count) = compositeSeries(:,i+1);
    count=count+1;
end

lineUpForMean = cat(3, LFP_Ex,LFP_In);
LFPs = mean(lineUpForMean, 3);
reg_name = ["S1","Thalamus","Insula","ACC","PFC","Amygdala","Hypothal", ...
    "LC","DRN","PB","NST"];
figure(2)

for i = 1:10
    subplot(5,2,i)
    plot(seriesSteps,LFPs(:,i+1))
    title(reg_name(:,i+1))
    xline(startStep)
    xline(endStep)
    xlim([leftEdge rightEdge])
    xticks([startStep, endStep])
    xticklabels([startSecs, endSecs])
    xlabel('time (s)')
end
f2 = gcf;
f2.Position = [300 300 1600 800];

% Cleaner version of S1 plot
figure(3)
hold on 
plot(seriesSteps, LFP)
ylim([-0.1 0.5])
xline(startStep)
xline(endStep)
xlim([leftEdge rightEdge])
xticks([startStep, endStep])
xticklabels([startSecs, endSecs])
ylim([bottomEdge topEdge])
xlabel('time (s)')
f3 = gcf;   
f3.Position = [400 400 2000 700];

% Package data for optional save, and create title for optional data and
% figure saves...
q.noiseval = noiseScaler;
q.w = p.w;
q.h = p.h;
q.a = p.a;
q.b = p.b;
q.epsilon = p.epsilon;
q.tau = p.tau;
q.dt = p.dt;
q.windowSecs = p.windowSecs;
q.threshold = p.threshold;
q.eucWeights = p.eucWeights;
q.noiseCHOICE = p.noiseCHOICE;
q.endtime = p.endtime;
q.n_runs = p.n_runs;
q.title = strcat(p.title, '_from_', string(startStep), '_for_', ...
    string(durStep), '_on_', string(Oix), '_', string(Tix), '_', ...
    string(Nix), '_ts');

% If required, save data...
if saveData
 save(q.title,'q')
end

% ...and plots
if savePlots
    saveas(f1,strcat(q.title, '.png'))
    saveas(f1,strcat(q.title, '.fig'))
    saveas(f1,strcat(q.title, '.eps'))
    saveas(f1,strcat(q.title, '.svg'))
    saveas(f2,strcat(q.title, '_others.png'))
    saveas(f2,strcat(q.title, '_others.fig'))    
    saveas(f2,strcat(q.title, '_others.eps'))    
    saveas(f2,strcat(q.title, '_others.svg'))    
    saveas(f3,strcat(q.title, '_clean.png'))
    saveas(f3,strcat(q.title, '_clean.fig'))
    saveas(f3,strcat(q.title, '_clean.eps'))
    saveas(f3,strcat(q.title, '_clean.svg'))
end