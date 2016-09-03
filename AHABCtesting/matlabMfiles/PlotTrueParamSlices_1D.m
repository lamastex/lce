% script to plot some 2-d marginal histograms
% 
% The input file gives a list of filenames that are used to make the
% graphics, one by one
% 
% Each image is saved separately
% 
% All images are named with the base name prepended and 
% saved in the current folder. 

clear functions
clear variables

%put the name of the input file here
filename = 'simdata2DGR_TrueParamNormalisedSliceSummaryStatMargs1D.txt';

%put the base name of the output files here
%this can include a path, in windows format, ie '\' path\file
outname = 'simdata2D';


%change/augment/shorten the sets of paramter and summary statistic names
% (these are used for labelling axes)
ssNames = {'heterozygosity''';'seg. sites'''};

f1 = @Function1DBoxesPlot;

[f d1] = ReadFilenamesAndSize(filename);



for i=1:size(d1)
    
    boxesFileName = f{1}{i};
    
    figure;

    h1 = gca;
    
    cla(h1);

    p = f1(boxesFileName, h1);
    set(get(h1,'XLabel'),'String',char(ssNames(d1(i),1)),'Interpreter', 'none');
   
    set(get(h1,'Title'),'String',boxesFileName,'Interpreter', 'none');
    
    
    outputfile = strcat(outname, 'TrueParamSlice', int2str(d1(i)), '.png');
    print ('-dpng', outputfile);
end

