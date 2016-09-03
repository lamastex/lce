% script to plot two 2d functions using real ranges
% 
% 

clear functions
clear variables

thesis = 1;

%put the name of the input file here
boxesFileName = 'simdata2DGR_Av_LeavesCherriesAvDepth_s_10_ns_25000_nr_10000_th_0.02.txt';
shapeFileName = 'simdata2DGR_AvPolygon_LeavesCherriesAvDepth_s_10_ns_25000_nr_10000_th_0.02.txt';

%parameter name
parameterName = 'theta''';
shapeParameterName = 'theta';
%parameterName = 'sigma';

%summary stat name
ssName = 'heterozygosity''';
shapeSSName = 'heterozygosity';

%put the base name of the output files here
%this can include a path, in windows format, ie '\' path\file
outname = 'simdata2D';

%change the view if necessary
v=[21 32]; 

%change the figure handle if necessary
f=figure;


f1 = @Function2DBoxesPlot;
f2 = @Function2DShapeBoxesPlot;


axisColour = [.3 .3 .3];
    
h1 = gca;
cla(h1);

p1 = f1(boxesFileName, h1);
set(get(h1,'XLabel'),'String',texlabel(parameterName));
set(get(h1,'YLabel'),'String',ssName,'Interpreter', 'none');
%set(h1,'View',[61 18]); % good for untransformed data
set(h1,'View',v); 
set(h1, 'FontSize',12, 'FontWeight','demi');
outputfile = strcat(outname, 'Average', '.png');
if (thesis == 0)
    set(get(h1,'Title'),'String',boxesFileName,'Interpreter', 'none');
else
    outputfile = strcat(outname, 'AverageThesis', '.png');
    set(h1, 'FontSize',12, 'FontWeight','demi');
    set(h1, 'XColor', axisColour, 'YColor', axisColour, 'ZColor', axisColour);
    set(get(h1,'XLabel'),'FontSize',12, 'FontWeight','demi');
    set(get(h1,'YLabel'),'FontSize',12, 'FontWeight','demi');
end


print ('-dpng', outputfile);


figure;

h2 = gca;
cla(h2);

p2 = f2(shapeFileName, h2);
set(get(h2,'XLabel'),'String',texlabel(shapeParameterName));
set(get(h2,'YLabel'),'String',shapeSSName,'Interpreter', 'none');
set(h2,'View',v); 

outputfile = strcat(outname, 'AveragePolygon', '.png');
if (thesis == 0)
    set(get(h2,'Title'),'String',shapeFileName,'Interpreter', 'none');
else
    outputfile = strcat(outname, 'AveragePolygonThesis', '.png');
    set(h2, 'FontSize',12, 'FontWeight','demi');
    set(h2, 'XColor', axisColour, 'YColor', axisColour, 'ZColor', axisColour);
    set(get(h2,'XLabel'),'FontSize',12, 'FontWeight','demi');
    set(get(h2,'YLabel'),'FontSize',12, 'FontWeight','demi');
end

print ('-dpng', outputfile);
