% script to plot some 2-d marginal histograms and slices at true params
% 
% The input file gives a list of filenames that are used to make the
% graphics, one by one
% 
% Each image is saved separately
% 
% All images are named with the base name prepended and 
% saved in the current folder. 

%put the name of the input file for the summary stat marginals here
filename = 'simdata4DGR_SummaryStatMargPolygons2D.txt';

%put the name of the input file for the true param slice summary stat marginals here
filenameTrueParams = 'simdata4DGR_TrueParamSliceSummaryStatMargPolygons2D.txt';

%put the base name of the output files here
%this can include a path, in windows format, ie '\' path\file
outname = 'simdata4DSmall';


%change the summary statistic names
% (these are used for labelling axes)
ssNames = {'heterozygosity';'seg. sites'};

%change the view if necessary
v = [-9 26];

clear functions
f1 = @Function2DShapeBoxesPlot;

[f d1 d2] = ReadFilenamesAndDims(filename);
[fTrue d1True d2True] = ReadFilenamesAndDims(filenameTrueParams);
assert(logical(size(d1,1) == size(d1True,1)));

for i=1:size(d1)
    
    assert(d1(i) == d1True(i));
    assert(d2(i) == d2True(i));

    boxesFileName = f{1}{i};
    boxesFileNameTrueParams = fTrue{1}{i};
    
    figure;

    h1 = gca;
    
    cla(h1);

    p = f1(boxesFileName, h1,[1,1,1],[.25,.25,.25],0,.25);
    set(get(h1,'XLabel'),'String',char(ssNames(d1(i),1)),'Interpreter', 'none');
    set(get(h1,'YLabel'),'String',char(ssNames(d2(i),1)),'Interpreter', 'none');
    
    set(h1,'View',v); 
   
    set(get(h1,'Title'),'String',...
        [{boxesFileName};{boxesFileNameTrueParams}],...
        'Interpreter', 'none');
    
    hold on
    p1 = f1(boxesFileNameTrueParams, h1);
    hold off
        
    outputfile = strcat(outname, 'SummaryStatMargAndTrueParamsMargSlicePolygon',...
                            int2str(d1(i)),'_', int2str(d2(i)), '.png');
    print ('-dpng', outputfile);
    
end

