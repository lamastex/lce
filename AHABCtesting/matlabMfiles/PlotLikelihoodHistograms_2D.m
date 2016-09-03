% script to plot some 2-d likelihood histograms
% 
% The input file gives a list of filenames that are used to make the
% graphics, one by one
% 
% Each image is saved separately and a gif of the images in sequence is 
% also made
% 
% All images and the gif are named with the base name prepended and 
% saved in the current folder.

%put the name of the input file here
filename = 'simdata4DGR_Likelihoods.txt';

%put the base name of the output files here
%this can include a path, in windows format, ie '\' path\file
outname = 'simdata4DSmall';

%change the parameter names
% these are used for labelling axes
% eg 'theta''' gives the effect of \theta^{\prime} in latex math
parameterNames = {'theta''','growth'''};

%change the gif name if necessary: it will automatically have the outname
% prepended using the format below
gifFilename = [outname 'Likelihood.gif'];

clear functions
f1 = @Function2DBoxesPlot;

[f s1] = ReadFilenamesAndSize(filename);

m = FunctionGetMaxHeight( f ,2);

figH = figure;

h1 = gca;

for i=1:size(s1)
    
    boxesFileName = f{1}{i};
    
    clf(figH)
    h1 = gca;
    set(h1,'NextPlot','replace');
   
    p = f1(boxesFileName, h1);
    
    set(get(h1,'XLabel'),'String',texlabel(parameterNames(1)));
    set(get(h1,'YLabel'),'String',texlabel(parameterNames(2)));
    set(get(h1,'ZLabel'),'String','Normalised Likelihood','Interpreter', 'none');

    set(h1,'ZLim',[0 m]);
    set(get(h1,'Title'),'String',boxesFileName,'Interpreter', 'none');
    
    set(h1,'View',[21 32]); 

    outputfile = strcat(outname, 'Likelihood', int2str(s1(i)), '.png');
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

