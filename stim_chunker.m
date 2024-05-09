% For fixed background input to NTS populations, and with other parameters 
% also fixed, this code will run a series of stochastic simulations, using 
% the same noise series, for different VNS amplitudes. There is a great
% deal of duplicated code from 'no_stim_chunker.m'
% 
% Combinations of different stimulation levels for both the excitatory and 
% inhibitory populations, although the example given here has all the 
% stimulation applied to the excitatory population of NTS. 
% 
% Due to RAM constraints, continuous runs are divided into epochs, each 
% starting with a consecutive noise seed. The output of this code is stored 
% in a compact format as arrays of tuples of varying length, recording the
% the starting timestep (for entire run) and duration of each seizure 
% event, as well as the initial conditions. From these we can calculate 
% seizure frequency and duration also reconstitute sections of the time 
% series of interest.
% 
% VNS is modelled as an idealised square wave, for which pulse width and 
% frequency can be specified. At any given point, VNS is calculated as a 
% function of the current time step for the whole run. This is added to the 
% NTS excitatory and inhibitory populations. To this end, the solver needs 
% to be given overall time step of the whole simulation.
% 
% The code can be parallelised (performing runs for more than parameter set 
% at a time) with MATLAB Parallel Toolbox. If not using this, change the 
% "parfor" loop, on line 226, to a "for" loop.

close all
clear

% can set number of parallel cores if using parallelisation
% parpool(44)

tic

%% SETUP PARAMETERS

maxExStim = 80;
nExStims = 76;

maxInStim = 0; 
nInStims = 1;
% (In this example, as in the paper, there is no VNS stimulation applied to
% the inhibitory population of the NTS)

% Background NTS(Ex) and (Inh) input values before any stimulation added.
baseEx = -0.9;
baseIn = -0.7;

% set number of epochs
nRuns = 100; % what was used for the examples in the paper.
% nRuns = 1; % to run quickly for testing

endtime = 100; % (simulated seconds for each section of epoch)  
% (NOTE cannot compare runs runs where sections are different lengths as 
% rng seeding will be different!! 100s used throughout paper)

% frequency of VNS pulse, in Hz (modelled as idealised square wave)
freq = 30; % default of 30 Hz used in paper.

% Pulsewidth of VNS pulse, in microseconds 
% (typically 250 or 500 in real-life clinical setting)
pulsewidth = 500; % default of 500 used in paper.

% add duty cycle information now for use if needed...
cycleOn = 300; % (seconds: default 30)
cycleOff = 300; % (default: 300)
% (Note: setting both to the same number results in no off periods. These
% are the parameters used in the paper)

%% Other important settings go here:

dt = 0.0001; % timestep: 0.0001s (100 microseconds) used throughout paper

% Select the dimensions contributing to Euclidean distance calculation
eucWeights = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; 
% (this setting considers only the S1 populations, as in the paper)

% width of window for moving average calculations
windowSecs = 2; % needs to be less than endtime
halfwindow = floor(windowSecs/(2*dt));
wyndoh = halfwindow * 2 +1; % easier if consistently an odd integer

% value used on moving average for seizure detection threshold
threshold = 0.15; 
% (0.15 suitable if using S1 populations to calculate euclidean distance. 
% Determined empirically. This setting was used in the paper)

% Calculate values of stim in Ex and Inh 
gapEx = maxExStim/(nExStims-1);
gapIn = maxInStim/(nInStims-1);
stimExValues = 0:gapEx:maxExStim;
stimInValues = 0:gapIn:maxInStim;

% get and modify standard parameters
p = read_default_params();

% (further information that may be passed between functions, and ultimately
% stored, can be put within the data structure, 'p')
p.startedAt = datetime;
p.dt = dt;
p.endtime = endtime;
p.nRuns = nRuns;

p.windowSecs = windowSecs;
p.threshold = threshold;
p.eucWeights = eucWeights;

p.freq = freq;
p.pw = pulsewidth;
p.cycleOff = cycleOff;
p.cycleOn = cycleOn;

% If needed, ammend default settings of any weights, offsets etc here:

% Templates for noiseINJECTED
noiseEX  = [1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];
noiseON  = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
noiseOFF = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
noiseNTS = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0];
noisePY  = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
noiseTC  = [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % (what we used)

% Noise settings...
nNoises = 1; % (sets the number of different noise levels to use)
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

%% Find precise FP and LC conditions 
% (in the absence of noise or stimulation). 

% The FP state is especially important for the detection of seizure-like 
% events. 

% We have considered doing array of FP/LC conditions for all parameter
% combinations in the sweep. This would have possibly made seizure
% detection more accurate, but is complicated by the problem that LC 
% may not exist at many of these combinations. Where it exists, however,
% the LC remains in a very similar position, as does the FP. For the
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

nearLC = zeros(1,22);
% (a starting state that reliably progresses to the limit cycle - when it
% exists - across a wide parameter range)

% Find equilibrated LC conditions once to save time later
[~,uLC] = vectorised_eulersolver ...
    (@(t,uLC)VNSfn_stoch_vec_Euler_stim(t,uLC,p), nearLC, dt, 10);
LCstart = uLC(end,:);

% reset for main loops
nSteps = int32(endtime/p.dt);
p.nSteps = nSteps;
p.nNoises = nNoises;

%% Set-up for Main Loop

% set overlaps (used for accurately caculating moving means)
origPreOverlap = uFP(end-halfwindow:end,:);

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

%% Main Loop

parfor lix = 1:listLength
% for lix = 1:listLength 
% (replace 'parfor' with 'for' if not parallelising)

    % state of system at end of run. Initiate at zero (no seizure)
    seizureFlag = 0;

    newInitState = exactFP; % (must be defined in parfor)
    preOverlap = origPreOverlap; % "
    localResult = zeros(2,0); % "
    tp = p; % (temporary version of broadcast variables)

    noiseScaler = sliceNoise(lix);
    tp.stimVal = [sliceOne(lix), sliceTwo(lix)];

    for seed = 1:nRuns
        
        % Work out how many timesteps have elapsed before current epoch
        tot_steps2epoch = (seed - 1) * nSteps;
    
        % Save initial state
        initStates(lix,seed,:) = newInitState;
        
        % Make some noise
        rng(seed, 'twister') 
        % (specifying algorithm keeps consistent between for and parfor
        % loops)
        baseNoise = randn(nSteps, 22);
        tp.noisevecs = baseNoise .* noiseCHOICE .* noiseScaler; 
        % (overwrites previous value)

        % additional info to pass to VNS_stim_fn 
        % (for for global timestep calculation)
        tp.epoch = seed;
  
        % send everything to the solver...
        [~,u] = vectorised_eulersolver(@(t,u)VNS_stim_fn(t,u,tp), ...
            newInitState, dt, endtime);

        % revise new initial states for next iteration
        newInitState = u(end,:);
       
        % calculate overlaps...
        rng(seed+1, 'twister')
        baseNoise = randn(nSteps, 22);
        tp.noisevecs = baseNoise .* noiseCHOICE .* noiseScaler; 
        % (overwrites previous value)
        tOverlap = halfwindow * tp.dt;
    
        % back to the solver...
        [~,postOverlap] = vectorised_eulersolver(@(t,u)VNS_stim_fn(t,u,tp), ...
            newInitState, dt, tOverlap);
        
        % Stitch together u its pre- and post- overlaps
        stitched = cat(1,preOverlap,u(2:end-1,:),postOverlap);
         
        % revise pre_overlap for next iteration
        preOverlap = u(end-halfwindow:end,:);
        
        % Analysis of segment (euclidean distances, moving averages, 
        % pruning, thresholding, coding):    

        % get euclidean distances...
        eucs = eucliser(stitched, eucWeights, exactFP);
        
        % get moving averages
        smooth_eucs = movmean(eucs,wyndoh);
        
        % trim ends to length...
        smooth_eucs = smooth_eucs(halfwindow+1:end-halfwindow-1,:);
    
        % apply thresholding
        isitoverthrshld = smooth_eucs > threshold;
    
        % A and B are two copies of 'isitoverthrshld', with the values of B
        % shifted by one position. 
        % Let C equal A minus B, then:
        %  a value of +1 implies a seizure start, 
        %  a value of -1 implies a seizure end, and
        %  a value of zero implies no change in seizure state. 
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
        fragmentResult = [zeros(2,0), [ons';lengths']]; 
        % (forces consistent dimensions)
    
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

    % assemble results here
    linearResults{lix} = localResult;

end % (of Main Loop)

%% Pack up and save results.

% Fold up linearResults...    
foldedResults(:) = linearResults(:);

% ...and initStates.
initStatesFolded(:) = initStates(:); 

p.tock = toc;

% Save the data tidily:

% tidy up bloat...
p = rmfield(p,'noisevecs');
p = rmfield(p,'stimVal');

% and add important data...
p.n_runs = nRuns;
p.endtime = endtime;
p.foldedResults = foldedResults;
p.init_states = initStatesFolded;
p.runUnder = mfilename;
p.finishedAt = datetime;
p.title = strcat("VNS_stim_output_", string(floor((now-738000)*1000)));
% (a general stem for for related plots etc. Therefore saved with data)

% save in appropriate subfolder...
save_dir = ['saved_output' filesep 'stim_chunker'];
[~,~] = mkdir (save_dir);
cd (save_dir)
save(p.title,'p')
cd(['..' filesep '..'])
