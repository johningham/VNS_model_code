This code accompanies the submitted paper: "Vagus nerve stimulation: Laying the groundwork for predictive network-based computer models" by John F. Ingham, Frances Hutchings, Yujiang Wang, Paolo Zuliani, and Peter N. Taylor. It was used in the development of the model and the preparation of the paper. The code was written by John F. Ingham and Frances Hutchings. 

The main directory contains six scripts that can each be run to produce different parts of figures used in the paper, and can be adapted to investigate the model further. The "Functions" directory contains functions called by these scripts. The "Data" directory contains data used in setting up the model, such as connection weights between polulations (in "VNSconnetivity.mat"), as well as examples of output produced by earlier stages of the pipeline to feed into later parts. 

By default, the scripts will set up a "saved_output" directory, within the main directory from which they are running, into which they will store their output to individualised sub-directories. This can obviously be changed from within the code. 

The six main scripts are as follows:

**composite_bifirc_fig_generator.m**

This code analyses the deterministic version of the model. It allows automated generation of sets of composite plots of bifurcation plots and time series plots for two-dimensional sweeps of the excitatory and inhibitory NTS background values, and to perform these for a series of altered values for any specified conection weight. Examples of this output were used in figure 2 of the our paper. Each subplot has a nominal pair of values for the fixed background input of the excitatory and inhibitory populations of the NTS. With each of these fixed in turn, the other is scanned through the ranges specified, producing a pair of bifurcation diagrams. Forward scans use red markers and backward scans use blue. When assembling these composite figures, many of the subplots will be reused, so the first part of the code produces these and stores them in a heirarchy of directories at a user specified location. 


**no_stim_chunker.m**

This code will run a series of stochastic simulations, for different combinations of constant input to NTS excitatory and inhibitory populations and noise levels. Other parameter combinations may be explored by changing "paramOne" and "paramTwo". For each combination of parameters, the same noise series are used.  

Due to RAM constraints, continuous runs are divided into epochs, each starting with a consecutive noise seed. The output of this code is stored in a compact format as arrays of tuples of varying length, recording the the starting timestep (for entire run) and duration of each seizure event, as well as the initial conditions. From this we can calculate seizure frequency and duration also reconstitute sections of the time series of interest. By default, the code is parallelised (performing runs for more than parameter set at a time) with MATLAB Parallel Toolbox. If not using this, change the "parfor" loop on line 201 to a "for" loop.


THIS CODE WILL TAKE HOURS TO RUN AS CURRENTLY SET UP (AS IN THE PAPER)!


**stim_chunker.m**

For fixed background input to NTS populations, and with other parameters also fixed, this code will run a series of stochastic simulations, using the same noise series, for different VNS amplitudes. There is a great deal of duplicated code from 'no_stim_chunker.m'. It is possible to apply combinations of different stimulation levels for both the excitatory and inhibitory populations, although the example given here has all the stimulation applied to the excitatory population of NTS, as was done for the examples in the paper.  

Due to RAM constraints, continuous runs are divided into epochs, each starting with a consecutive noise seed. The output of this code is stored in a compact format as arrays of tuples of varying length, recording the the starting timestep (for entire run) and duration of each seizure event, as well as the initial conditions. From these we can calculate seizure frequency and duration also reconstitute sections of the time series of interest. 

VNS is modelled as an idealised square wave, for which pulse width and frequency can be specified. At any given point, VNS is calculated as a function of the current time step for the whole run. This is added to the NTS excitatory and inhibitory populations. To this end, the solver needs to be given overall time step of the whole simulation. By default, the code is parallelised (performing runs for more than parameter set at a time) with MATLAB Parallel Toolbox. If not using this, change the "parfor" loop, on line 229, to a "for" loop.

THIS CODE WILL TAKE HOURS TO RUN AS CURRENTLY SET UP (AS IN THE PAPER)!


**stim_reconstitutor.m**

Takes an output file from 'stim_chunker.m' and reconstitutes parts of the time series that may be of interest for a specified start point, duration, and combination of parameters. Time series without VNS can be obtained from 'stim_chunker.m', where stimulation amplitude is set to zero. The code outputs not only time series of the mean S1 populations' time series, but also the euclidean distance from the FP, and the smoothed version of this used in the seizure detector. This code was used on similar data to produce the supplementary figure that explained the seizure detection process. The code also outputs a figure of time series for the remaining regions over the same period with the same parameters; a version of the S1 time series (without superimposed euclidean distance information); and a final figure showing a phase plot of the S1 populations, again similar to the one used in the supplementary figure.


**one_D_no_stim_plotter.m**

Takes output file from 'no_stim_chunker.m' and produces plots of both seizure frequency, and proportion of total time in spent in seizure, against the background input to the excitatory NTS population. The latter plot was used in Figure 3 of the paper.


**one_D_stim_plotter.m**

Takes output file from 'stim_chunker.m' and produces plots of both seizure frequency, and proportion of total time in spent in seizure, against the amplitude of the VNS applied to the the excitatory NTS population (ultimately, the only one we applied it to). The latter plot was used in Figure 4 of the paper.


# Requirements:
MATLAB (by Mathworks). Runs on versions 2022b through 2024a

Parallel Computing Toolbox (by Mathworks)

Signal Processing Toolbox (by Mathworks)

Digital Signal Processing Using MATLAB, Version 1.0.0.0 (by John Proakis)

