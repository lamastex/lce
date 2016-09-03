% script to plot leaves and log posteriors
% 
%

%put the name of the input file here
boxesFileName = 'LeafLogPost1simdata2DGR_s_10_ns_25000_nr_6000000_th_0.02.log';
%change the figure handle if necessary
figure(10);

clear functions
f1 = @FunctionLeafAndLogPosteriorPlot;

h1 = gca;
cla(h1);

p = f1(boxesFileName, h1);
%set(p(1),'XLim',[150000 300000]);
%set(p(2),'XLim',[150000 300000]);


