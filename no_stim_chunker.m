% This code will run a series of stochastic simulations, for different 
% combinations of constant input to NTS excitatory and inhibitory 
% populations and noise levels. Other parameters may be scanned by changing 
% "paramOne" and "paramTwo" For each combination of parameters, the same 
% noise series are used. 
% 
% Due to RAM constraints, continuous runs are divided into epochs, each 
% starting with a consecutive noise seed. The output of this code is stored 
% in a compact format as the starting timestep and duration of each seizure 
% event, as well as the initial conditions. From this we can calculate 
% seizure frequency and duration also reconstitute sections of the time 
% series of interest.
% 
% The code can be parallelised (performing runs for more than parameter set 
% at a time) with MATLAB Parallel Toolbox. If not using this, change the 
% "parfor" loop on line *** to a "for" loop.

close all
clear

% can set number of parallel cores if using parallelisation
% parpool(4)

tic

%% SETUP PARAMETERS

paramOne = 21; % a single value for p.h(i) or pair ([i,j]) for p.w(i,j)
% (usually the NTS(Py) offset -- p.h(21) and its values)
paramOneValues = -1:0.01:0; % what was used in the paper

paramTwo = 22; % a single value for p.h(i) or pair ([i,j]) for p.w(i,j)
% (usually the NTS(Inh) offset -- p.h(22) and its values)
paramTwoValues = 0; % specify single value when not wishing to sweep this parameter

nParamOne = length(paramOneValues);
nParamTwo = length(paramTwoValues);

% set number of epochs
nRuns = 100; % what was used for the examples in the paper.
% nRuns = 1; % to run quickly for testing

endtime = 100; % (simulated seconds for each section of epoch)  
% (NOTE cannot compare runs runs where sections are different lengths as 
% rng seeding will be different!! 100s used throughout paper)

%% Other important settings go here:

dt = 0.0001; % timestep: 0.0001s (100 microseconds) used throughout paper

% Select the dimensions contributing to Euclidean distance calculation
eucWeights = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; 
% (this setting considers only the S1 populations, as in the paper)

% width of window for moving average calculations
windowSecs = 2; % needs to be less than endtime (not much sense if not anyway)
halfwindow = floor(windowSecs/(2*dt));
wyndoh = halfwindow * 2 + 1; % easier if consistently an odd integer

% value used on moving average for seizure detection threshold
threshold = 0.15; 
% (0.15 suitable if using S1 populations to calculate euclidean distance. 
% Determined empirically. This setting was used in the paper)

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

p.paramOne = paramOne;
p.ParamOneValues = paramOneValues;
p.nParamOne = nParamOne;
p.paramTwo = paramTwo;
p.paramTwoValues = paramTwoValues;
p.nParamTwo = nParamTwo;

% If needed, ammend defuault settings of any weights, offsets etc here:

% Templates for noiseINJECTED
noiseEX  = [1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];
noiseON  = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
noiseOFF = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
noiseNTS = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0];
noisePY  = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
noiseTC  = [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % (what we used)

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

%% Find precise FP and LC conditions 
% (in the absence of noise or stimulation). 

% The LC parameter is especially important for the detection of
% seizure_like events. 

% We have considered doing array of FP/LC conditions all parameter
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
listLength = nParamOne * nParamTwo * nNoises;
listParamTuples = zeros(3,listLength);
for nix = 1:nNoises
    for tix = 1:nParamTwo
        for oix = 1:nParamOne
            tuple = [paramOneValues(oix), paramTwoValues(tix), noiseScalers(nix)];
            lix = (nix-1)*nParamTwo*nParamOne + (tix-1)*nParamOne + oix;
            listParamTuples(:,lix) = tuple;
        end
    end
end

% Set up cell array to store results for each parameter set. This has the
% same linear indexing as the tuples above.
linearResults = repmat({int64(zeros(2,0))}, 1, listLength);
% initialise results array 
foldedResults = repmat({int64(zeros(2,0))},nParamOne, nParamTwo, nNoises);

% create array to store initial states of all runs (for ease of recreating
% particular time segments). 
initStates = zeros(listLength, nRuns, 22); 
% Generated with all parameter combinations in one dimension for ease of 
% use in Parfor loop, but folded up later for ease of use.make this up 
% empty to fill later...
initStatesFolded = zeros(nParamOne, nParamTwo, nNoises, nRuns, 22);

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
    
    % deal with both offsets and weights as parameters...
    if isscalar(paramOne) % i.e. an offset value
        tp.h(paramOne) = sliceOne(lix); 
    else % it should be a weight in the form [x,y]
        tp.w(paramOne(1), paramOne(2)) = sliceOne(lix);
    end
    if isscalar(paramTwo) % i.e. an offset value
        tp.h(paramTwo) = sliceTwo(lix); 
    else % it should be a weight in the form [x,y]
        tp.w(paramTwo(1), paramTwo(2)) = sliceTwo(lix);
    end

    noiseScaler = sliceNoise(lix);

    for seed = 1:nRuns
        
        % work out how many timesteps have elapsed before current epoch
        tot_steps2epoch = (seed - 1) * nSteps;
    
        % save initial state
        initStates(lix,seed,:) = newInitState;
        
        % Make some noise
        rng(seed, 'twister') 
        % (specifying algorithm keeps consistent between for and parfor
        % loops)
        baseNoise = randn(nSteps, 22);
        tp.noisevecs = baseNoise .* noiseCHOICE .* noiseScaler; 
        % (overwrites previous value)
 
        % send everything to the solver...
        [~,u] = vectorised_eulersolver(@(t,u) ...
            VNSfn_stoch_vec_Euler_stim(t,u,tp), newInitState, dt, endtime);

        % revise new initial states for next iteration
        newInitState = u(end,:);
       
        % calculate overlaps...
        rng(seed+1,'twister')
        baseNoise = randn(nSteps, 22);
        tp.noisevecs = baseNoise .* noiseCHOICE .* noiseScaler; % (overwrites previous val)
        tOverlap = halfwindow * tp.dt;
    
        % back to the solver...
        [~,postOverlap] = vectorised_eulersolver(@(t,u) ...
            VNSfn_stoch_vec_Euler_stim(t,u,tp), newInitState, dt, tOverlap);
        
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

    % results assembled here
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

p.title = strcat("VNS_no_stim_output_", string(floor((now-738000)*1000)));
% (a general stem for for related plots etc. Therefore saved with data)
p.runUnder = mfilename;
p.finishedAt = datetime;

save(p.title,'p')