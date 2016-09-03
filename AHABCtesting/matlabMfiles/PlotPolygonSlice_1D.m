% script to plot one 2d shape using real ranges
% plotting the 2d box and real range as a 3-d box
% 

%put the name of the input file here
boxesFileName = 'simdata2DGR_LikelihoodPolygon_s_10_ns_25000_nr_10000_th_0.02_20.txt';
%change the figure handle if necessary
figure;

clear functions
f1 = @Function1DShapeBoxesPlot;

h1 = gca;
cla(h1);

p = f1(boxesFileName, h1);
color_edge = [0.65 0.65 0.65];    % edge colour
set(p,'EdgeColor', color_edge);
%set(h1,'Xlim',[0.5 1]);
%set(h1,'Ylim',[0.5 1]);
set(get(h1,'Title'),'String',boxesFileName,'Interpreter', 'none');
set(get(h1,'XLabel'),'String',texlabel('theta'));


