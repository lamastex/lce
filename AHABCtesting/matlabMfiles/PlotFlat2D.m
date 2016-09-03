% script to plot one 2d function using real ranges
% 
% 

%put the name of the input file here
%boxesFileName = 'mixtureCarvedGRAverage.txt';
%boxesFileName = 'mixtureCarvedSingleAverage.txt';
%boxesFileName = 'mixtureTransformedCarvedGRAverage.txt';
%boxesFileName = 'mixtureTransformedBestOne.txt';
boxesFileName = 'simdata2DGR_Av_s_10_ns_25000_nr_10000_th_0.02.txt';

%change the figure handle if necessary
figure;

clear functions
f1 = @Function2DFlatBoxesPlot;

h1 = gca;
cla(h1);

p = f1(boxesFileName, h1);
set(get(h1,'XLabel'),'String',texlabel('theta'));
set(get(h1,'YLabel'),'String','h','Interpreter', 'none');
%set(h1,'View',[61 18]); % good for untransformed data

set(get(h1,'Title'),'String',boxesFileName,'Interpreter', 'none');
