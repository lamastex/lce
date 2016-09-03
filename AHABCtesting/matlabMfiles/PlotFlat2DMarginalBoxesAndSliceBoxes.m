% script to plot some 2-d marginal histograms
% 
% The input file gives a list of filenames that are used to make the
% graphics, one by one
% 
% Each image is saved separately
% 
% All images are named with the base name prepended and 
% saved in the current folder. 

%put the name of the input file here
filename = 'simdata4DGR_Margs.txt';

%put the names of the file(s) to use as slice boxes here
%assumes one slice on each marginal histogram listed in the input filename
sliceBoxesFileNames = {'sliceBoxes988572554_0.0574807.txt',...
                        'sliceBoxes3573308719_0.0873993.txt',...
                        'sliceBoxes2183992837_0.0574807.txt',...
                        'sliceBoxes3403369844_0.0873993.txt'};


%put the base name of the output files here
%this can include a path, in windows format, ie '\' path\file
outname = 'simdata4DSmall';


%change/augment/shorten the sets of paramter and summary statistic names
% (these are used for labelling axes)
parameterNames = {'theta''';'growth'''};
ssNames = {'heterozygosity''';'seg. sites'''};

%change the view if necessary
v=[90 90];

clear functions
f1 = @Function2DFlatBoxesPlot;
f2 = @Function2DFlatSliceBoxesPlot;

[f d1 d2] = ReadFilenamesAndDims(filename);

for i=1:size(d1)
    
    boxesFileName = f{1}{i};
    sliceBoxesFileName = sliceBoxesFileNames{i};
    
    slDiv1 = regexp(sliceBoxesFileName, '_', 'split');
    slDiv2 = regexp(slDiv1{2}, '\.txt', 'split');
    slpt = slDiv2{1};
   
    figure;

    h1 = gca;
    
    cla(h1);

    p = f1(boxesFileName, h1, [1 1 1], [.8 .8 .8], 0, .2);
    
    hold on
    p2 = f2(sliceBoxesFileName, h1);
    hold off
    
    set(get(h1,'XLabel'),'String',texlabel(char(parameterNames(d1(i),1))));
    set(get(h1,'YLabel'),'String',char(ssNames(d2(i),1)),'Interpreter', 'none');
    %set(h1,'View',[61 18]); % good for untransformed data
    set(h1,'View',v); 

    gtitle=[{boxesFileName};{sliceBoxesFileName}];
    set(get(h1,'Title'),'String',gtitle,'Interpreter', 'none');
    
    slpt = regexprep(slpt, '\.', 'p')
    outputfile = strcat(outname, 'FlatMargSlice_', int2str(d1(i)),...
                '_', int2str(d2(i)),...
                '_', slpt, '.png');
    print ('-dpng', outputfile);
end

