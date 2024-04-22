% For fixed background input to NTS populations, and with other parameters 
% also fixed, this code will run  a series of stochastic simulations, using 
% the same noise seeds, for different VNS amplitudes. Effective stimulation 
% levels can be specified for both the excitatory and inhibitory
% populations, and all combinations tested, although the example given here
% has all the stimulation applied to the excitatory population of NTS.
%
% VNS is modelled as an idealised square wave, for which pulse width and 
% frequency can be specified. At any given point, VNS is calcuated as a 
% function of the current time step for the whole run. This is added to the
% NTS excitatory and inhibitory populations. To this end, the solver needs 
% to keep an eye on the overall time step of the whole simulation.
%
% The output of this code is stored in a compact format that allows
% calculation of seizure frequency and duration and also allows for
% infrequent events to be readily reconstitiuted 
 
close all
clear

% can set number of parallel cores if using parallelisation
% parpool(44)

tic

% SETUP PARAMS HERE!!!*****************************************************
maxExStim = 80;
nExStims = 76;

maxInStim = 0; 
nInStims = 1;
% (In this example, as in the paper, there is no VNS stimulation applied to
% the inhibitory population of the NTS)

% Background NTS(Ex) and (Inh) input values before stim added
baseEx = -0.9;
baseIn = -0.7;

% set number of epochs
nRuns = 100;

endtime = 100; % (simulated seconds for each section of epoch)  
% (NOTE cannot compare runs runs where sections are different lengths as rng
% seeding will be different!!)

% frequency of VNS pulse in Hz (modelled as idealised square wave)
freq = 30; % default of 30 Hz, as used in paper.

% pulsewidth of VNS pulse, in microseconds 
% (250 or 500 in real-life clinical setting)
pulsewidth = 500; % as used in paper.

% add duty cycle information now for use if needed...
cycleOn = 300; % (seconds: default 30)
cycleOff = 300; %(default: 300)
% (Note: setting both to the same number results in no off periods. This
% are the parameters used in the paper)

% value used on moving average for seizure detection threshold
threshold = 0.15; 
% (0.15 suitable if just sticking with S1 to calculate euclidean distance. 
% Determined empirically. This setting was used in the paper)

% *************************************************************************


% Other important settings go here:

dt = 0.0001; % (timestep)

% use to switch dimensions contributing to Euclidean distance
% calculation...
% eucWeights = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0]; % (everything but NTS)
eucWeights = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % (just S1)

% width of window for moving average calculations
windowSecs = 2; % needs to be less than endtime
halfwindow = floor(windowSecs/(2*dt));
wyndoh = halfwindow*2 +1; % easier if consistently an odd integer

% Calculate values of stim in Ex and Inh 
gapEx = maxExStim/(nExStims-1);
gapIn = maxInStim/(nInStims-1);
stimExValues = 0:gapEx:maxExStim;
stimInValues = 0:gapIn:maxInStim;

% get and modify standard parameters
p = read_default_params();

% (anything else to pass to the function via 'p')
p.dt = dt;
p.windowSecs = windowSecs;
p.threshold = threshold;
p.eucWeights = eucWeights;
p.endtime = endtime;
p.freq = freq;
p.pw = pulsewidth;
p.cycleOff = cycleOff;
p.cycleOn = cycleOn;

% (then ammend any weights, offsets, noise scaling etc)
p.w(4,21) = 0.08; % default = 0.08, but connectivity.mat value is 0.01
p.w(4,15) = 0.01; % default = 0.01, but connectivity.mat value is ZERO
p.h(21) = -1.5; 
% (default = -0.9 but for now this is set here in order to calculate a
% consistent value for the FP, based on earlier work)

% templates for noiseINJECTED
noiseEX = repmat([1,0],1,11); 
noiseON = ones(1,22);
noiseOFF = zeros(1,22);
noiseNTS = noiseOFF; noiseNTS(21) = 1;
noisePY = noiseOFF; noisePY(1) = 1;
noiseCustom1 = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
noiseCustom2 = [1 0 0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];
noiseCustom3 = [1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

% Noise settings...
nNoises = 1; % (sets the number of different noise levels that will be used)
noiseMaxScaler = 0.72;
noiseMinScaler = 0;
% This is set to produce the single noise scaler (effectively the standard
% deviation of the noise) of 0.72 that is used in the paper.
noiseCHOICE = noisePY;

% allow for single value of noise scaler (equal to Max)
if nNoises == 1
    noiseScalers = noiseMaxScaler;
else
    noiseInterval = (noiseMaxScaler - noiseMinScaler)/ (nNoises - 1);
    noiseScalers = noiseMinScaler:noiseInterval:noiseMaxScaler;
end

p.noiseCHOICE = noiseCHOICE;
p.noiseScalers = noiseScalers;

% Initial condition sets:
% (S1_PY, SI_IN, TC, RE, INS_EX, INS_IN, ACC_EX, ACC_IN, PFC_EX, PFC_IN, 
% Amy_Ex, Amy_In, Hyp_Ex, Hyp_In, LoC_Ex, LoC_In, DRN_Ex, DRN_In, PB_Ex,
% PB_In, STN_Ex, STN_In)

nearFP = [0.1724,0.1787,-0.0818,0.2775,0.0724,0.0787,0.0724,0.0787,...
    0.0724,0.0787,0.0724,0.0787,0.0724,0.0787,0.0724,0.0787,0.0724,...
    0.0787,0.0724,0.0787,0.1724,0.1787];
   
nearLC = zeros(1,22);

%% Find precise FP conditions (for current params, no noise, no stim) for measurements later

% (We have debated doing array of FP/LC conditions all parameter
% combinations in the sweep. This would have possibly made seizure
% detection more accurate, but is complicated by the problem that LC or FP
% will not exist at many of these combinations. A such, we are keeping the
% original, lower dimension functions for this purpose. A bit inelegant,
% but pragmatic...)
nSteps = int32(10/p.dt + 1);  % (runs for 10 secs, done once in single dimension, so low overhead)
p.stimVal = [0,0]; % Leave it!
p.noisevecs = zeros(nSteps, 22);
[~,uFP] = vectorised_eulersolver(@(t,uFP)VNSfn_stoch_vec_Euler_stim(t,uFP,p), nearFP, dt, 10);
exactFP = uFP(end,:);

% Find equilibrated LC conditions (for current params, no noise, no stim) once to save time later
[~,uLC] = vectorised_eulersolver(@(t,uLC)VNSfn_stoch_vec_Euler_stim(t,uLC,p), nearLC, dt, 10);
LCstart = uLC(end,:);

% set overlaps (used for accurately caculating moving means)
origPreOverlap = uFP(end-halfwindow:end,:);

% (reset n_steps for main loops)
nSteps = int32(endtime/p.dt);
p.nSteps = nSteps;
p.nNoises = nNoises;

% Set up the list of tuples
listLength = nExStims * nInStims * nNoises;
listParamTuples = zeros(3,listLength);
for nix = 1:nNoises
    for iix = 1:nInStims
        for eix = 1:nExStims
            tuple = [stimExValues(eix), stimInValues(iix), noiseScalers(nix)];
            lix = (nix-1)*nInStims*nExStims + (iix-1)*nExStims + eix;
            listParamTuples(:,lix) = tuple;
        end
    end
end

% Set up cell array to store results for each parameter set. This has the
% same linear indexing as the tuples above.
linearResults = repmat({int64(zeros(2,0))}, 1, listLength);
% initialise results array  
foldedResults = repmat({int64(zeros(2,0))},nExStims, nInStims, nNoises);

% create array to store initial states of all runs (for ease of recreating
% particular time segments). 
initStates = zeros(listLength, nRuns, 22); 
% Generated with all parameter combinations in one dimension for ease of 
% use in Parfor loop, but folded up later for ease of use.make this up 
% empty to fill later...
initStatesFolded = zeros(nExStims, nInStims, nNoises, nRuns, 22);

sliceOne = listParamTuples(1,:);
sliceTwo = listParamTuples(2,:);
sliceNoise = listParamTuples(3,:);

%% PARFOR LOOP STARTS HERE!
% parfor lix = 1:listLength
for lix = 1:listLength % replace with FOR if not parallelising

    % state of system at end of run. Initiate at zero (no seizure)
    seizureFlag = 0;

    newInitState = exactFP; % (must be defined in parfor)
    preOverlap = origPreOverlap; % "
    localResult = zeros(2,0); % "
    tp = p; % (temporary version of broadcast variables)

    noiseScaler = sliceNoise(lix);
    tp.stimVal = [sliceOne(lix), sliceTwo(lix)];

    for seed = 1:nRuns
        
        % work out how many timesteps have elapsed before current epoch
        tot_steps2epoch = (seed - 1) * nSteps;
    
        % save initial state
        initStates(lix,seed,:) = newInitState;
        
        % make some noise
        rng(seed, 'twister') % specifying algorithm keeps consistent between for and parfor loops
        baseNoise = randn(nSteps, 22);
        tp.noisevecs = baseNoise .* noiseCHOICE .* noiseScaler; % (overwrites previous val)

        % additional info to pass to function
        tp.epoch = seed;
    
        % send everything to the solver...
        [~,u] = vectorised_eulersolver(@(t,u)VNS_stim_fn(t,u,tp), newInitState, dt, endtime);

        % revise new initial states for next iteration
        newInitState = u(end,:);
       
        % calculate overlaps...
        rng(seed+1, 'twister')
        baseNoise = randn(nSteps, 22);
        tp.noisevecs = baseNoise .* noiseCHOICE .* noiseScaler; % (overwrites previous val)
        tOverlap = halfwindow * tp.dt;
    
        % back to the solver...
        [~,postOverlap] = vectorised_eulersolver(@(t,u)VNS_stim_fn(t,u,tp), ...
            newInitState, dt, tOverlap);
        
        % Stitch together u its pre- and post- overlaps
        stitched = cat(1,preOverlap,u(2:end-1,:),postOverlap);
         
        % revise pre_overlap for next iteration
        preOverlap = u(end-halfwindow:end,:);
        
        % Analysis of segment (euclidean distances, moving averages, pruning, thresholding,
        % coding):    
        % get euclidean distances...
        eucs = eucliser(stitched, eucWeights, exactFP);
        
        % get moving averages
        % (currently moving mean but could convolve with hamming/hann window etc?)
        smooth_eucs = movmean(eucs,wyndoh);
        
        % trim ends to length...
        % (would need to cut off more if using convolution)
        smooth_eucs = smooth_eucs(halfwindow+1:end-halfwindow-1,:);
    
        % apply thresholding
        isitoverthrshld = smooth_eucs > threshold;
    
        % (New method...)
        A = isitoverthrshld;
        B = A;
        B(2:end) = B(1:end-1);
        B(1) = 0;
        C = int8(A - B);
        
        ons = int32(find(C==1));
        offs = int32(find(C==-1));
        
        if A(end)==1
            offs(end+1,1)=length(A)+1; % (always end final run of ones)
        end
        
        lengths = offs - ons;
        fragmentResult = [zeros(2,0), [ons';lengths']]; % (forces consistent dims)
    
        % (join contiguous seizure elements)
        if seizureFlag == 1 && isitoverthrshld(1) == 1 
            localResult(2,end) = localResult(2,end) + lengths(1);
            fragmentResult = fragmentResult(:,2:end);
        end

        % convert from epoch to global timestep...        
        fragmentResult(1,:) = fragmentResult(1,:) + tot_steps2epoch;

        localResult = cat(2, localResult, fragmentResult);
    
        % update seizure flag            
        seizureFlag = isitoverthrshld(end);

    end % (of seed/epoch)

    % assemble results here!!
    linearResults{lix} = localResult;

end % (of parfor loop!)

% Fold up linearResults...    
foldedResults(:) = linearResults(:);

% and initStates...
initStatesFolded(:) = initStates(:); 

p.tock = toc;
toc

% Save the data tidily:

% tidy up bloat...
p = rmfield(p,'noisevecs');
p = rmfield(p,'stimVal');

% and add important data...
p.n_runs = nRuns;
p.endtime = endtime;
p.foldedResults = foldedResults;
p.init_states = initStatesFolded;

p.title = strcat("VNS_stim_output_", string(floor((now-738000)*1000)));
% (p.title is a general stem for for related plots etc. Therefore saved with data)
fulltitle = [p.title, '.mat'];

save(p.title,'p')
