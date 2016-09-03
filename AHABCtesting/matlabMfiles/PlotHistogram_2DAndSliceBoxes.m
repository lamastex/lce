% script to plot one 2d function using real ranges
% 
% 

%put the name of the input file here
boxesFileName = 'simdata2DGR_Av_s_10_ns_25000_nr_10000_th_0.02.txt';
sliceBoxesFileNames = {'sliceBoxes1596879076_0.527764.txt',...
                        'sliceBoxes3150612281_0.13718.txt',...
                        'sliceBoxes3274023777_0.058268.txt',...
                        'sliceBoxes676046140_-0.0710424.txt',...
                        'sliceBoxes675538646_0.113971.txt'};
%put the base name of the output files here
%this can include a path, in windows format, ie '\' path\file
outname = 'simdata2DGR';

%change the parameter name
% this is used for labelling the graphic (in latex math form)
% use for example 'theta''' to get theta^{prime} effect
parameterName = 'theta''';

%summary stat name
% this is used for labelling the graphic
% use for example 'h''' to get h^{prime} effect
ssName = 'h''';

%change the figure handle if necessary
figure;

clear functions
f1 = @Function2DBoxesPlot;
f2 = @Function2DSliceBoxesPlot;

h1 = gca;
cla(h1);

p1 = f1(boxesFileName, h1, [1 1 1], [.6 .6 .6], 0, 1);
set(get(h1,'XLabel'),'String',texlabel(parameterName));
set(get(h1,'YLabel'),'String',ssName,'Interpreter', 'none');
set(h1,'View',[21 32]); 

gtitle={boxesFileName};

p = []; % array of patch graphic objects for the slices

slpts = '';

hold on;
for i=1:size(sliceBoxesFileNames,2)
    
    sliceBoxesFileName = sliceBoxesFileNames{i};
    
    slDiv1 = regexp(sliceBoxesFileName, '_', 'split');
    slDiv2 = regexp(slDiv1{2}, '\.txt', 'split');
    slpts = [slpts, '_', slDiv2{1}];
       
    p2 = f2(sliceBoxesFileName, h1);
    p=[p p2];

    gtitle = [gtitle;{sliceBoxesFileName}];
end
set(get(h1,'Title'),'String',gtitle,'Interpreter', 'none');

hold off;

slpts = regexprep(slpts, '\.', 'p');
outputfile = strcat(outname, 'HistSlices', slpts, '.png');
    print ('-dpng', outputfile);