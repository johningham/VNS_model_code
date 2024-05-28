% Takes an output file from 'stim_chunker.m' and reconstitutes parts of the
% time series that may be of interest for a specified start point, 
% duration, and combination of parameters. 
% 
% The code outputs not only time series of the mean S1 populations' time
% series, but also the euclidean distance from the FP, and the smoothed
% version of this used in the seizure detector. This code was used on
% similar data to produce the supplementary figure that explained the 
% seizure detection process.

% The code also outputs a figure of time series for the remaining regions
% over the same period with the same parameters, a version of the time
% time raw S1 time series (without euclidean distance information), and a 
% final figure showing a phase plot of the S1 populations, again similar to
% the one used in the supplementary figure.

close all
clear

% ensure that we have the working directory matches the location of script
main_folder = fileparts(which(mfilename));
cd(main_folder)

% add all subdirectories to the path
addpath(genpath(main_folder))

load('VNS_stim_output_1147915.mat') % data used in final paper

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
f1 = figure(1);
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
% f1 = gcf;
f1.Position = [200 200 2000 700];

%% Plotting other regions
f2 = figure(2);

stepsToPlot = size(seriesSteps,2);
ex_pops = zeros(stepsToPlot, 11);
in_pops = zeros(stepsToPlot, 11);

for i = 1:10
    ex_pops(:,i) = compositeSeries(:,2*i + 1);
    in_pops(:,i) = compositeSeries(:,2*i + 2);
end

lineUpForMean = cat(3, ex_pops,in_pops);
LFPs = mean(lineUpForMean, 3);
reg_name = ["Thalamus","Insula","ACC","PFC","Amygdala", "Hypothalamus", ...
    "LC","DRN","PB","NST"];

for i = 1:10
    subplot(5,2,i)
    plot(seriesSteps,LFPs(:,i))
    title(reg_name(:,i))
    xline(startStep)
    xline(endStep)
    xlim([leftEdge rightEdge])
    xticks([startStep, endStep])
    xticklabels([startSecs, endSecs])
    xlabel('time (s)')
end

% f2 = gcf;
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

%% Phase plot

start = 1;
fin = 100000;

f4 = figure(4);
hold on

h1 = plot(compositeSeries(start:fin,1), compositeSeries(start:fin,2));
h1.Color = [.4 .4 .4];
h1.LineWidth = 1;

% Show Fixed Point as blue dot
h2 = plot(exactFP(1), exactFP(2), 'o');
h2.Color = 'b';

h2.MarkerSize = 8;
h2.MarkerFaceColor = 'b';

ylabel('S1 Py')
xlabel('S1 Inh')

% Show Limit Cycle in with red dashed line
load('shortCleanSeries.mat')
% (presaved segment of LC under our standard conditions, with no noise and
% no VNS)

howFarRound = 3370;
% (ensures dashed line not overdrawn by using only one period of LC)

h3 = plot(clean_composite_series(1:howFarRound,1), ...
    clean_composite_series(1:howFarRound,2));
h3.Color = [1,0,0,0.5];
h3.LineWidth = 3;
h3.LineStyle = ":";

legend('','FP','LC', Location = 'northwest')


%% If required, save data...
if saveData
    save_dir = ['saved_output' filesep 'reconstitutor' filesep 'data'];
    [~,~] = mkdir (save_dir);
    cd (save_dir)
    save(q.title,'q')
    cd(['..' filesep '..' filesep '..'])
end

% ...and plots
if savePlots
    save_dir = ['saved_output' filesep 'reconstitutor' filesep 'plots'];
    [~,~] = mkdir (save_dir);
    cd (save_dir)
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
    saveas(f4,strcat(q.title, '_phase.png'))
    saveas(f4,strcat(q.title, '_phase.fig'))
    saveas(f4,strcat(q.title, '_phase.eps'))
    saveas(f4,strcat(q.title, '_phase.svg'))
    cd(['..' filesep '..' filesep '..'])
end
