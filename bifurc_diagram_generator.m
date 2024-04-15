% This code analyses the deterministic version of the model. It allows
% automated generation of sets of composite plots for two-dimensional 
% sweeps of the excitatory and inhibitory NTS background values, and to
% perform these for a series of values for different conection weights.
%
% As set up, the code takes some considerable time to generate a set of
% composite figures for all cominations of 21 excitatory and 21 inhibitory
% NTS values, a total of 441. It performs these for only a single value of
% the NTS>TC connection weight, but more could be specified. The run time
% for this task is in the order of

close all
clear
% set location for saved output
loc = '/Users/john/Documents/VNS_output/composite_bifurc_figs';
%% ************************TWEAK VARIABLES HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% The example settings will 
wt2vary = [4, 21]; % e.g.[4,21] for NTS>TC
blurb = 'NTS2RE';
% ('blurb' prefixes autogenerated title. It should refer to the connection
% weight being varied. It also ensures that it has its own directory into
% which it will store the output figures)

weights = 0.01; 
% (lists or ranges may be used to scan over a number of different values of
% the desired connection weight)

% Below is an example of an extensive range of parameter values
ex_min = -5; 
ex_max = 5; 
ex_step = 0.5; 

inh_min = -8; 
inh_max = 2; 
inh_step = 0.5; 


%% ************************************************************************

% Get standard parameters and tweak as needed...
p = read_default_params();

p.h(3) = -2; % change from default for bistable ODE45
base_params = p; % save copy to reset from as needed.

% (Autocalculate steps for outer parameters)
exes = ex_min:ex_step:ex_max;
inhs = inh_min:inh_step:inh_max;

% (Auto calculate for inidvidual bifurcation sweeps)
ex_range = ex_max:-0.01:ex_min; 
inh_range = inh_max:-0.01:inh_min;

% Ensure the existence of the folders we will need.
fldr = [loc, blurb];
mkdir(fldr) % top level. composite plots here as well as subdirectories
cd(fldr)
mkdir mat % contained in <fldr>. composite plots saved here as matlab figs
mkdir parts % also in <fldr>. individual plots that make the composite plots
cd parts 
mkdir mat % contained in <fldr>\parts\. matlab fig versions of individual plots

for wix = 1:length(weights) % (everything goes inside this loop!)

    % Set strength of relevant connection 
    connection_wt = weights(wix);
    p.w(wt2vary(1),wt2vary(2)) = connection_wt; 
        
    % Loop and create NTSEx sweeps for various NTSInh. Save to 'parts'.
    param_to_change = 21; % (NTS excitatory)
    paramrange = ex_range; 
    for iix = 1:length(inhs)
        NTS_inh = inhs(iix);
        p.h(22) = NTS_inh;
        fig = Bifurcation_VNS_takes_params(p, param_to_change, paramrange, 2);
        ylabel('Min, Max S1') 
        xlabel('NTS Ex input')
        fig.XLim = [ex_min, (ex_max + (ex_max - ex_min)*0.1)]; % (add a bit of space on the right)
        filename = ['ExSweep', blurb, '=', num2str(connection_wt), ' NTS(inh)=', num2str(NTS_inh)];
        title(filename) % (just for its individual existence)
    
        % (needs saving to parts and parts\mat!)
        oldFolder = cd(fldr);
        cd parts
        saveas(fig,[filename, '.png'])
        cd mat
        saveas(fig,[filename, '.fig'])
        cd(oldFolder) 
        close all
    end
    p.h(22) = base_params.h(22);

    % Loop and create NTS_Inh sweeps for various NTS_Ex. Save to 'parts'
    param_to_change = 22; % (NTS inhibitory)
    paramrange = inh_range; 
    for iix = 1:length(exes)
        NTS_ex = exes(iix);
        p.h(21) = NTS_ex;
        fig = Bifurcation_VNS_takes_params(p, param_to_change, paramrange, 2);
        ylabel('Min, Max S1') 
        xlabel('NTS Inh input')
        fig.XLim = [inh_min, (inh_max + (inh_max - inh_min)*0.1) + 0.0001]; 
        % (add a bit of space on the right, and avoids error if min+max)
        filename = ['InhSweep', blurb, '=', num2str(connection_wt), ' NTS(ex)=', num2str(NTS_ex)];
        title(filename) % (just for its individual existance)
    
        % % (needs saving to parts and parts\mat!)
        % oldFolder = cd(fldr);
        % cd parts
        % saveas(fig,[filename, '.png'])
        % cd mat
        % saveas(fig,[filename, '.fig'])
        % cd(oldFolder)
        % close all
    end
    p.h(21) = base_params.h(21);
    
     
    %% The nested for loops within which time series are generated, and
    % everything is assembled, named and saved
    
    for eix = 1:length(exes)
        for iix = 1:length(inhs)

            NTS_ex = exes(eix);
            NTS_inh = inhs(iix); 
    
            close all

            % fetch ax1 and ax2 from file. Alter labels and titles as needed.
            % Subplot each. 
            
            filename1 = ['ExSweep', blurb, '=', num2str(connection_wt),...
                ' NTS(inh)=', num2str(NTS_inh), '.fig'];
            filename2 = ['InhSweep', blurb, '=', num2str(connection_wt),...
                ' NTS(ex)=', num2str(NTS_ex) '.fig'];

            % oldFolder = cd(fldr);
            % cd parts
            % cd mat
            % 
            % openfig(filename1);
            % figure(1)
            % ax1 = subplot(2, 2, 1, gca);
            % title('')
            % xline(NTS_ex, 'b--', {'current value    '})
            % 
            % openfig(filename2);
            % figure(2)
            % ax2 = subplot(2, 2, 2, gca);
            % title('')
            % ylabel('')
            % xline(NTS_inh, 'b--', {'current value    '})
            % 
            % cd(oldFolder)
            

            % DO CORESPONDING TIME SERIES PLOTS.
            p.h(21) = NTS_ex;
            p.h(22) = NTS_inh;
            % Get figure for time series at FP (if present) by calling plot_time_series_fn.
            % (initial condition near the fixed point)
            init_cond = [0.1724,0.1787,-0.0818,0.2775,0.0724,0.0787,...
                0.0724,0.0787,0.0724,0.0787,0.0724,0.0787,0.0724,0.0787,...
                0.0724,0.0787,0.0724,0.0787,0.0724,0.0787,0.1724,0.1787];
            
            ax3 = subplot(2,2,3,plot_time_series_fn(p, init_cond));
                xlabel('Time (sec)','FontSize',20)
                ylabel('S1','FontSize',20)
                title('"FP"')
                
            % Get figure for time series at LC (if present) by calling plot_time_series_fn.
            % (quick and dirty way to start close to seizure)
            init_cond = zeros(1,22); 
            
            ax4 = subplot(2,2,4,plot_time_series_fn(p, init_cond));
                xlabel('Time (sec)','FontSize',20)
                ylabel('')
                title('"LC"')
            
            % Get all Y-axes on same scale

            h1 = findobj(ax1, 'Type', 'scatter');
            yd1 = get(h1, 'YData');
            hc1 = horzcat(yd1{:});
            ax1_min = min(hc1);
            ax1_max = max(hc1);

            h2 = findobj(ax2, 'Type', 'scatter');
            yd2 = get(h2, 'YData');
            hc2 = horzcat(yd2{:});
            ax2_min = min(hc2);
            ax2_max = max(hc2);

            ax1_lim = [ax1_min, ax1_max];
            ax2_lim = [ax2_min, ax2_max];
           
            ylims = vertcat(ax1_lim, ax2_lim);%, ax3.YLim, ax4.YLim); 
            % (For now ignore YLims of time series)
            ylimsMin = min(ylims(:,1));
            ylimsMax = max(ylims(:,2));
            NuLim = [ylimsMin, (ylimsMax + 0.1)]; % (add some headroom for title)
            ax1.YLim = NuLim;
            ax2.YLim = NuLim;
            ax3.YLim = NuLim;
            ax4.YLim = NuLim;
            
            % Put together composite figure
            composite_figure = figure;
            titl = [blurb,' = ', num2str(connection_wt),'  NTS(Ex) = ',...
                num2str(NTS_ex), '  NTS(Inh) = ', num2str(NTS_inh)];
            sgtitle (titl)
            copyobj(ax1, composite_figure)
            copyobj(ax2, composite_figure)
            copyobj(ax3, composite_figure)
            copyobj(ax4, composite_figure)
    
            %Size figure consistently (?needed)
%             fig.Position = [1000 200 1200 1100];
    
            % Add code to save current figure in relevant place with relevant title.
            % For now... 
            oldFolder = cd(fldr);
            saveas(composite_figure,[titl, '.png'])
            cd mat
            saveas(composite_figure,[titl, '.fig'])
            cd(oldFolder)
            
            close all
        
        end
    end
end

