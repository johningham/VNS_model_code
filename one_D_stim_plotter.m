% Takes output file from 


clear
close all

% load output file from which 
load('VNS_stim_output_1147915.mat') %(as used in final paper)



dotColour = [0.6350 0.0780 0.1840];

% Specify which noise scaler value to use (if more than one)
noiseIx = 1;
    
tot_timesteps = p.n_runs * p.nSteps; 
ex_vals = p.stimExVals;
in_vals = p.stimInVals;
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
        pc_dur_seizures(ex,in) = (sum(cell_mat(2,:)) * 100) ./ cast(tot_timesteps,'double');
    end
end

seizuresPerHour = n_seizures ./ totHours;

f1 = figure(1);
s1 = plot(ex_vals,seizuresPerHour,'o');
xlabel('Excitatory VNS component')
ylabel('Number of seizures per hour')
s1.Color = dotColour;
s1.MarkerSize = 8;
s1.MarkerFaceColor = dotColour;
ylim([0, 1])
f1.Position = [700 400 900 200];
% saveas(gcf,strcat(p.title, 'Numb.png'))
% saveas(gcf,strcat(p.title, 'Numb.fig'))
% saveas(gcf,strcat(p.title, 'Numb.eps'))

f2 = figure(2);
s2 = plot(ex_vals,pc_dur_seizures,'o');
xlabel('Excitatory VNS component')
ylabel('% of total time spent in seizure')
s2.Color = dotColour;
s2.MarkerSize = 8;
s2.MarkerFaceColor = dotColour;
ylim([0, 0.4])
f2.Position = [700 100 900 200];
% saveas(gcf,strcat(p.title, 'Pc.png'))
% saveas(gcf,strcat(p.title, 'Pc.fig'))
% saveas(gcf,strcat(p.title, 'Pc.eps'))
