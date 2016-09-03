% script to plot some 1-d likelihood polygon shapes with coverage regions
% as well
% 
% The input file gives a list of filenames that are used to make the
% graphics, one by one
% 
% Each image is saved separately and a gif of the images in sequence is 
% also made
% 
% All images and the gif are named with the base name prepended and 
% saved in the current folder.

clear functions
clear variables

%put the name of the likelihood input file here
filenameLikelihoods = 'simdata2DGR_LikelihoodPolygons.txt';

%put the name of the coverage input file here
filenameCoverages = 'simdata2DGR_CoveragePolygons.txt';

%put the name of the coverage area input file here
filenameCoverageAreas = 'simdata2DGR_CoverageAreaPolygons.txt';

%put the base name of the output files here
%this can include a path, in windows format, ie '\' path\file
outname = 'SimData2D';

%change the parameter name
% this is used for labelling axes
% eg 'theta' gives the effect of \theta in latex math
parameterName = 'theta';

%change the gif name if necessary: it will automatically have the outname
% prepended using the format below
gifFilename = [outname 'LikelihoodWitCoveragePolygon.gif'];

f1 = @Function1DShapeBoxesPlot;
f2 = @Function1DShapeCoverageBoxesPlot;

[f s1] = ReadFilenamesAndSize(filenameLikelihoods);

[fCov s1Cov] = ReadFilenamesAndSize(filenameCoverages);

[fCovArea s1CovArea] = ReadFilenamesAndSize(filenameCoverages);

m = FunctionGetMaxHeight( f,1 );

figH = figure;

h1 = gca;

for i=1:size(s1)
    
    boxesFileName = f{1}{i};
    covFileName = fCov{1}{i};
    
    clf(figH)
    h1 = gca;
    %set(h1,'NextPlot','replace');
   
    p1 = f1(boxesFileName, h1);
    hold on
    p2 = f2(covFileName, h1, [0 0 1]);
    hold off
    set(get(h1,'XLabel'),'String',texlabel(parameterName));
    set(get(h1,'YLabel'),'String','Normalised Likelihood','Interpreter', 'none');

    set(h1,'YLim',[0 m]);
    
    set(get(h1,'Title'),'String',boxesFileName,'Interpreter', 'none');
    
    outputfile = strcat(outname, 'LikelihoodPolygonWithCoverage', int2str(s1(i)), '.png');
    print ('-dpng', outputfile);
    
    set(get(h1,'Title'),'String','');
    
    frame = getframe(gcf);
   
   im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1;
    imwrite(imind,cm,gifFilename,'gif', 'Loopcount',1);
    else
    imwrite(imind,cm,gifFilename,'gif','WriteMode','append');
    end
end

