% Takes an output file from stim_chunker.m and reconstitutes part of the
% time series for a specified interval and combination of parameters.

close all
clear

load('VNS_stim_output_1147915.mat') % data for use in final paper

savePlots = false;
saveData = false;

% choose from available param combinations and specify start time and
% duration
Oix = 2;
Tix = 1;
Nix = 1;
startStep =  245628404;
durStep = 433520;

startmargin = 50000;
endmargin = 100000;

topEdge = 0.5;
bottomEdge = -0.1;

% Start needs to be no bigger than pre-calc timesteps 

endStep = startStep + durStep;

% unpack p
threshold = p.threshold;
dt = p.dt;
endtime = p.endtime;

stimExVal = p.stimExVals(Oix);
stimInVal = p.stimInVals(Tix);
noiseval = p.noiseScalers(Nix);

windsecs = p.windowSecs;
halfwindow = floor(windsecs/(2*dt));
wyndoh = halfwindow*2 +1; % (easier if consistently an odd integer)
clear("nStps")
nStps = cast(p.nSteps,'double');

leftEdge = startStep - startmargin; % (x-limits of plot)
rightEdge = endStep+endmargin;

verystart = cast((startStep - startmargin - halfwindow),'double');
veryend = cast((startStep + durStep + endmargin + halfwindow), 'double');

vsepoch = cast(floor(verystart/nStps),'double') + 1; % epochs numbered from one. Seed = epoch 
vslocaltimestep = verystart - ((vsepoch-1) * nStps); 

veepoch = cast((floor(veryend/nStps) + 1), 'double'); 
velocaltimestep = veryend - ((veepoch-1) * nStps); 

initsOTNER = p.init_states; % (paramOne, paramTwo, Noise, Epoch, Region)
initsER = permute((initsOTNER(Oix,Tix,Nix,:,:)),[4 5 1 2 3]);


% BORROWED CODE (cut, pasted and hacked about)
%% Find precise FP conditions (for current params, no noise, no stim) for measurements later

% (We have debated doing array of FP/LC conditions all parameter
% combinations in the sweep. This would have possibly made seizure
% detection more accurate, but is complicated by the problem that LC or FP
% will not exist at many of these combinations. A such, we are keeping the
% original, lower dimension functions for this purpose. A bit inelegant,
% but pragmatic...)

% Initial condition sets:
% (S1_PY, SI_IN, TC, RE, INS_EX, INS_IN, ACC_EX, ACC_IN, PFC_EX, PFC_IN, 
% Amy_Ex, Amy_In, Hyp_Ex, Hyp_In, LC_Ex, LC_In, DRN_Ex, DRN_In, PB_Ex,
% PB_In, STN_Ex, STN_In)

% Find exact FP (runs for 10 secs, done once in single dimension, so low overhead)
nearFP = [0.1724,0.1787,-0.0818,0.2775,0.0724,0.0787,0.0724,0.0787,...
    0.0724,0.0787,0.0724,0.0787,0.0724,0.0787,0.0724,0.0787,0.0724,...
    0.0787,0.0724,0.0787,0.1724,0.1787];
nStps = int32(10/p.dt + 1);  
p.stimVal = [0,0]; % Leave it!
p.noisevecs = zeros(nStps, 22);
[~,uFP] = vectorised_eulersolver(@(t,uFP)VNSfn_stoch_vec_Euler_stim(t,uFP,p), nearFP, dt, 10);
exactFP = uFP(end,:);

% Find equilibrated LC conditions (for current params, no noise, no stim) once to save time later
nearLC = zeros(1,22);
[~,uLC] = vectorised_eulersolver(@(t,uLC)VNSfn_stoch_vec_Euler_stim(t,uLC,p), nearLC, dt, 10);
LCstart = uLC(end,:);

% noise 
p.AffIn = p.noiseCHOICE * noiseval;

% initiate long time series
compositeSeries = zeros(0,22);

% handle calculations that need to start before time point zero
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
    p.noisevecs = randn(nStps, 22) .* p.AffIn; % (overwrites previous val)
    
    % get the relevant initial conditions to start the epoch
    thisinit = initsER((seed),:);

    % send to solver
    [~,u] = vectorised_eulersolver(@(t,u)VNS_stim_fn(t,u,p), thisinit, ... 
        dt, endtime-dt);

    % append result to long time series
    compositeSeries = cat(1,compositeSeries,u); 
end

% Calculate a little further on for safety...
seed = seed + 1;
p.epoch = seed;
rng(seed, 'twister')
p.noisevecs = randn(nStps, 22) .* p.AffIn; % (overwrites previous val)

% get the relevant initial conditions to start the epoch
thisinit = compositeSeries(end,:);

% send to solver
[~,u] = vectorised_eulersolver(@(t,u)VNSfn_stoch_vec_Euler_stim(t,u,p), ...
    thisinit, dt, endtime-dt);

% append result to long time series
compositeSeries = cat(1,compositeSeries,u); 

% Calculate LFP (mean of S1 populations)
LFP = mean(compositeSeries(:,1:2),2);

% calculate euclidean distance (using freestanding function)
eucs = eucliser(compositeSeries, p.eucWeights, exactFP);

% get moving averages
smooth_eucs = movmean(eucs,wyndoh);

% trim ends to length...
% smooth_eucs = smooth_eucs(halfwindow+1:end-halfwindow-1,:);

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
% showing all regions, there are 11, but we are plotting the cortical S1
% traces in figure 1 so missing that.
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
q.noiseval = noiseval;
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

% If required,save data...
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
