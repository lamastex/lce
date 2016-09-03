% script to plot one 2d function in top down view, ie just using domain
% boxes, and add coloured 'slice' boxes to this
% 
% 

clear functions
clear variables

%put the name of the input files here
boxesFileName = 'simdata2DGR_Av_s_10_ns_25000_nr_600000_th_0.02.txt';
sliceBoxesFileNames = {'simdata2DGR_SliceBoxes8_-0.748834.txt',...
                        'simdata2DGR_SliceBoxes1_0.504256.txt',...
                        'simdata2DGR_SliceBoxes7_-0.0388621.txt',...
                        'simdata2DGR_SliceBoxes2_0.106375.txt',...
                        'simdata2DGR_SliceBoxes10_-0.209093.txt'};

%put the base name of the output files here
%this can include a path, in windows format, ie '\' path\file
outname = 'simdata2DGR';

%parameter name
parameterName = 'theta''';

%summary stat name
ssName = 'h''';


%change the figure handle if necessary
figure;

f1 = @Function2DFlatBoxesPlot;
f2 = @Function2DFlatSliceBoxesPlot;

h1 = gca;
cla(h1);

p1 = f1(boxesFileName, h1,[1 1 1], [.6 .6 .6], 0, 1);
set(h1,'View',[0 90]); 
set(get(h1,'XLabel'),'String',texlabel(parameterName));
set(get(h1,'YLabel'),'String',ssName,'Interpreter', 'none');

gtitle={boxesFileName};
slpts = '';

p = []; % array of patch graphic objects for the slices

hold on;
for i=1:size(sliceBoxesFileNames,2)
    
    sliceBoxesFileName = sliceBoxesFileNames{i};
    
    slDiv1 = regexp(sliceBoxesFileName, '_', 'split');
    slDiv2 = regexp(slDiv1{2}, 'SliceBoxes', 'split');
    slpts = [slpts, '_', slDiv2{2}];
       
    p2 = f2(sliceBoxesFileName, h1);
    p=[p p2];

    gtitle = [gtitle;{sliceBoxesFileName}];
end
set(get(h1,'Title'),'String',gtitle,'Interpreter', 'none');

hold off;
slpts = regexprep(slpts, '\.', 'p');
outputfile = strcat(outname, 'FlatBoxSlices', slpts, '.png');
    print ('-dpng', outputfile);
