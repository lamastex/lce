% script to plot one 2d function using real ranges
% 
% 

%put the name of the input file here
boxesFileName = 'normalsdGR_Av_nr_10000_sd_1.0.txt';

%parameter name
%parameterName = 'theta''';
parameterName = 'sigma';

%summary stat name
ssName = 'x';

%change the figure handle if necessary
figure;

clear functions
f1 = @Function2DBoxesPlot;

h1 = gca;
cla(h1);

p = f1(boxesFileName, h1);
set(get(h1,'XLabel'),'String',texlabel(parameterName));
set(get(h1,'YLabel'),'String',ssName,'Interpreter', 'none');
%set(h1,'View',[61 18]); % good for untransformed data
set(h1,'View',[21 32]); 

set(get(h1,'Title'),'String',boxesFileName,'Interpreter', 'none');
