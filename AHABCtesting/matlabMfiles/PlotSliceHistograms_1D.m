% script to plot some 1-d slice histograms
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

%put the name of the input file here
filename = 'simdata2DGR_Slices.txt';

%put the base name of the output files here
%this can include a path, in windows format, ie '\' path\file
outname = 'simdata2DGR';

%change the parameter name
% this is used for labelling axes
% eg 'theta''' gives the effect of \theta^{\prime} in latex math
parameterName = 'theta''';

%change the gif name if necessary: it will automatically have the outname
% prepended using the format below
gifFilename = [outname 'Slices.gif'];

f1 = @Function1DBoxesPlot;

[f v1] = ReadFilenamesAndVal(filename);

m = FunctionGetMaxHeight( f ,2);

figH = figure;

h1 = gca;


for i=1:size(v1)
    
    boxesFileName = f{1}{i};
    
    clf(figH)
    h1 = gca;
    set(h1,'NextPlot','replace');
   
    p = f1(boxesFileName, h1);
     
    set(get(h1,'XLabel'),'String',texlabel(parameterName));
    set(h1,'YLim',[0 m]);

    ystring = ['Un-normalised conditional density at ', num2str(v1(i),'% 10.2f')];
    set(get(h1,'YLabel'),'String',ystring,'Interpreter', 'none');
    set(get(h1,'Title'),'String',boxesFileName,'Interpreter', 'none');
    outputfile = strcat(outname, 'Slice', int2str(i), '.png');
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

