% script to plot leaves and log posteriors
% 
%

%put the name of the input file here
boxesFileName = 'LeafLogPost1simdata2DGR_s_10_ns_25000_nr_6000000_th_0.02.log';
%change the figure handle if necessary
figure;

clear functions
f1 = @FunctionLeafAndLogPosterior2ChainsPlot;

h1 = gca;
cla(h1);

p = f1(boxesFileName, h1);
set(get(h1,'Title'),'String',boxesFileName);

