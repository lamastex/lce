function MyTajDPlots(UncondCol, CondCol, CondAlpha, PlotTitle,FileName)
pvs = dlmread(FileName, '\t');

% 'number of unique x^*s seen in data is'
length(unique(pvs(:,2)))

IndxPosS1 = find(pvs(:,4)>0);


pvS1 = pvs(IndxPosS1,:);

SizeIndxPosS1 = size(pvS1)

%UnCondAlpha = qthSampleQuantile(0.05, sort(pvS1(:,UncondCol)))
%CondAlpha = qthSampleQuantile(0.05, sort(pvS1(:,CondCol)))

clf
%figure
colormap(hot)
%colormap(gray)
%loglog([1e-3,1],[1e-3,1])
subplot(1,2,1)
title(strcat('',PlotTitle))

hold on
Alpha=0.05
plot([Alpha,Alpha],[Alpha,1])
plot([Alpha,Alpha],[0.0001,Alpha])
plot([Alpha,0.001],[Alpha,Alpha])
plot([Alpha,1],[Alpha,Alpha])
scatter(pvS1(:,UncondCol),pvS1(:,CondCol),8,pvS1(:,3),'o','filled')
xlabel('Unconditional p-values'); ylabel('Conditional p-values');colorbar
%one-sided test??
%scatter(pvS1(:,9),pvS1(:,10),12,pvS1(:,3),'o','filled')
hold off
% how many are rejected by un/conditional test


IndxPosSRejUnCond = find(pvs(:,4)>0 & pvs(:,UncondCol) < Alpha);
NumPosRejUnCond = length(IndxPosSRejUnCond)



IndxPosSRejCond = find(pvs(:,4)>0 & pvs(:,CondCol) < CondAlpha);
NumPosSRejCond = length(IndxPosSRejCond)

IndxPosSRejBoth = find(pvs(:,4)>0 & pvs(:,UncondCol) < Alpha & pvs(:,CondCol) < CondAlpha);
NumPosSRejBoth = length(IndxPosSRejBoth)

FalseRejects = find(pvs(:,4)>0 & pvs(:,UncondCol) < Alpha & pvs(:,CondCol) >= CondAlpha);

NumberFalseRejects = length(FalseRejects)

%sfsFileName = '../work/Data/SFS_n_10.t_100.r_0.g_0';
%sfsFileName = '../../sim/MS_Simul_Data/n_10.t_100.g.r_0_10_100.m_100000/SFS_n_10.t_100.g_0.r_0.m_100000.reps_10000';
%sfs = dlmread(sfsFileName, ' ');
%size(sfs)
%sfs(FalseRejects,:)

% EDF of p-values and 5% Confidence Band
%figure
subplot(1,2,2)
title(strcat('',PlotTitle))
hold on;

Alpha=0.05; % set alpha to 5% for instance

% Get the x and y coordinates of SampleSize-based ECDF in x1 and y1 and
plot([-0.0001,1],[0,1])
% plot the first ECDF using the function ECDF
[x1 y1] = ECDF(pvS1(:,UncondCol),0,0.01,0.01);
stairs(x1,y1,'b');
SampleSize=size(pvS1(:,UncondCol),1)
Epsn = sqrt((1/(2*SampleSize))*log(2/Alpha)); % epsilon_n for the confidence band
stairs(x1,max(y1-Epsn,zeros(1,length(y1))),'g'); % lower band plot
stairs(x1,min(y1+Epsn,ones(1,length(y1))),'g'); % upper band plot
axis([0 1 0 1]);
%axis square;

% plot the first ECDF using the function ECDF
[x1 y1] = ECDF(pvS1(:,CondCol),0,0.01,0.01);
stairs(x1,y1,'r');
SampleSize=size(pvS1(:,CondCol),1);
Epsn = sqrt((1/(2*SampleSize))*log(2/Alpha)); % epsilon_n for the confidence band
stairs(x1,max(y1-Epsn,zeros(1,length(y1))),'g'); % lower band plot
stairs(x1,min(y1+Epsn,ones(1,length(y1))),'g'); % upper band plot
xlabel('Uncond. or Cond. p-values'); ylabel('EMF');
%LabelString=['n=' num2str(SampleSize)];
%text(0.75,0.05,LabelString)

MeanS = mean(pvs(:,4))
hold off;