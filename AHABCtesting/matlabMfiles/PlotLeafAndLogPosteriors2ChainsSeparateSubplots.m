% script to plot leaves and log posteriors on 2 separate subplots
% 
%

%put the name of the input file here
%boxesFileName = 'LeafLogPost0AutoGRScalarsKeepBest.log';
boxesFileName = 'LeafLogPostTrace_Leaves_Logpost0AutoGRScalarsGuassian6DKeepApartTwoFlags.log';

%put the base name of the output files here
%this can include a path, in windows format, ie '\' path\file
outname = 'Gaussian6D';

%change the figure handle if necessary
figure;

clear functions
f1 = @FunctionLeaf2ChainsPlot;
f2 = @FunctionLogPosterior2ChainsPlot;

h = gca;
cla(h);

h1 = subplot(1,2,1,'replace');
p1 = f1(boxesFileName, h1);
h2 = subplot(1,2,2,'replace');
p2 = f2(boxesFileName, h2);
set(get(h1,'Title'),'String',boxesFileName);

outputfile = strcat(outname, 'LeafAndLogPostTrace', '.png')
    
print ('-dpng', outputfile);