% Takes output file from 'no_stim_chunker.m' and produces plots of both
% seizure frequency, and proportion of total time in spent in seizure,
% against the background input to the excitatory NTS population. The latter 
% plot was used in Figure 3 of the paper.

clear
close all

% ensure that we have the working directory matches the location of script
main_folder = fileparts(which(mfilename));
cd(main_folder)

% add all subdirectories to the path
addpath(genpath(main_folder))

% set "true" to save plots
savePlots = false;

% Specify index of the noise scaler value to use (only one in this dataset)
noiseIx = 1;

% load output file from which to draw plots... 
load('VNS_no_stim_output_1146671.mat') % (as used in final paper)

dotColour = [0.9290 0.6940 0.1250]; % (a kind of mustard yellow?)

tot_timesteps = p.n_runs * p.nSteps; 
ex_vals = p.ParamOneValues;
in_vals = p.paramTwoValues;
n_exes = length (ex_vals);
n_inhs = length (in_vals);
n_seizures = zeros(n_exes,n_inhs);
pc_dur_seizures = zeros(n_exes,n_inhs);

totSecs = p.n_runs * p.endtime;
totHours = totSecs/3600;

for ex = 1:n_exes
    for in = 1:n_inhs
        cell_mat = p.foldedResults{ex,in,noiseIx};
        n_seizures(ex,in) = length(cell_mat);
        pc_dur_seizures(ex,in) = (sum(cell_mat(2,:)) * 100) ./ ...
            cast(tot_timesteps,'double');
    end
end

seizuresPerHour = n_seizures ./ totHours;

f1 = figure(1);
s1 = plot(ex_vals,seizuresPerHour,'o');
xlabel('NTS_{PY} input')
ylabel('Number of seizures per hour')
s1.Color = dotColour;
s1.MarkerSize = 8;
s1.MarkerFaceColor = dotColour;
ylim([0, 1])
f1.Position = [700 400 900 200];

f2 = figure(2);
s2 = plot(ex_vals,pc_dur_seizures,'o');
xlabel('NTS_{PY} input')
ylabel('% of total time spent in seizure')
s2.Color = dotColour;
s2.MarkerSize = 8;
s2.MarkerFaceColor = dotColour;
ylim([0, 0.4])
f2.Position = [700 100 900 200];

% Save plots if flag set to true
if savePlots
    save_dir = ['saved_output' filesep 'no_stim_plotter'];
    [~,~] = mkdir (save_dir);
    cd (save_dir)
    saveas(f1,strcat(p.title, 'Numb.png'))
    saveas(f1,strcat(p.title, 'Numb.fig'))
    saveas(f1,strcat(p.title, 'Numb.eps'))
    saveas(f2,strcat(p.title, 'Pc.png'))
    saveas(f2,strcat(p.title, 'Pc.fig'))
    saveas(f2,strcat(p.title, 'Pc.eps'))
    cd(['..' filesep '..'])
end
