% script to plot one 2d function using real ranges
% 
% 

%put the name of the input file here
%boxesFileName = 'mixtureCarvedGRSlice.txt';
%boxesFileName = 'mixtureBestOneMCMCSlice.txt.txt';
%boxesFileName = 'mixtureTransformedBestOneSlice.txt';
boxesFileName = 'simdata2DGR_Likelihood_s_10_ns_25000_nr_10000_th_0.02_10.txt';
%change the figure handle if necessary
figure;

clear functions
f1 = @Function1DBoxesPlot;

h1 = gca;
cla(h1);

p = f1(boxesFileName, h1);
set(get(h1,'XLabel'),'String',texlabel('theta'));
set(get(h1,'Title'),'String',boxesFileName,'Interpreter', 'none');
