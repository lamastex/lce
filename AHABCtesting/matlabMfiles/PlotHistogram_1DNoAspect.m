% script to plot one 2d function using real ranges
% 
% 

%put the name of the input file here
boxesFileName = 'normalsdGR_Likelihood_nr_100000_sd_1.0.txt';

%change the figure handle if necessary
figure;

clear functions
f1 = @Function1DBoxesPlot;

h1 = gca;
cla(h1);

p = f1(boxesFileName, h1);
set(get(h1,'XLabel'),'String',texlabel('sigma'));

set(get(h1,'Title'),'String',boxesFileName,'Interpreter', 'none');
