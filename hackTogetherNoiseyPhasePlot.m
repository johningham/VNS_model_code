close all
clear
load('shortNoisySeries.mat')

start = 40000;
fin = 110000;

f = figure(4);
hold on

h1 = plot(compositeSeries(start:fin,1), compositeSeries(start:fin,2));
h1.Color = [.4 .4 .4];
h1.LineWidth = 1;

% Show Fixed Point as blue dot
h2 = plot(exactFP(1), exactFP(2), 'o');
h2.Color = 'b';

h2.MarkerSize = 8;
h2.MarkerFaceColor = 'b';

% legend('LC','FP', Location = 'north')

ylabel('S1 Py')
xlabel('S1 Inh')

% Show Limit Cycle in with red dashed line
load('shortCleanSeries.mat') 
% (presaved segment of LC under our standard conditions, with no noise and
% no VNS)

howFarRound = 3370;
% (ensures dashed line not overdrawn by using only one period of LC)

h3 = plot(compositeSeries(1:howFarRound,1), compositeSeries(1:howFarRound,2));
h3.Color = [1,0,0,0.5];
h3.LineWidth = 3;
h3.LineStyle = ":";