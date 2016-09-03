FileName = '../work/Tests/SFS_n_10.t_1_10_100.r_0.g_0_mcE4.pvalues.txt';

pvs = dlmread(FileName, '\t');

% 'number of unique x^*s seen in data is'
length(unique(pvs(:,2)))
%scatter(pvs(:,7),pvs(:,8), 'k')

IndxPosS1 = find(pvs(:,4)>0 & mod(pvs(:,1),3)==1);
IndxPosS10 = find(pvs(:,4)>0 & mod(pvs(:,1),3)==2);
IndxPosS100 = find(pvs(:,4)>0 & mod(pvs(:,1),3)==0);


pvS1 = pvs(IndxPosS1,:);
pvS10 = pvs(IndxPosS10,:);
pvS100 = pvs(IndxPosS100,:);

size(pvS1)
size(pvS10)
size(pvS100)

clf
%figure
%colormap(hot)
colormap(gray)
p1=loglog([1e-4,1],[1e-4,1])
axis([0.0001 1 0.0001 1])
hold on
Alpha=0.05
p2=plot([Alpha,Alpha],[Alpha,1])
p3=plot([Alpha,Alpha],[0.0001,Alpha])
p4=plot([Alpha,0.0001],[Alpha,Alpha])
p5=plot([Alpha,1],[Alpha,Alpha])
p6=scatter(pvS1(:,7),pvS1(:,8),10,pvS1(:,3),'o')
p7=scatter(pvS10(:,7),pvS10(:,8),10,pvS10(:,3),'+')
p8=scatter(pvS100(:,7),pvS100(:,8),10,pvS100(:,3),'*')
legend([p6,p7,p8],'theta=1','theta=10','theta=100')
title('Unconditional versus conditional p-values of Tajima`s D test of Standard Neutrality')
xlabel('Unconditional p-values'); ylabel('Conditional p-values');colorbar
IndxPosSRejUnCond1 = find(pvs(:,4)>0 & pvs(:,7) < Alpha & mod(pvs(:,1),3)==1);
IndxPosSRejCond1 = find(pvs(:,4)>0 & pvs(:,8) < Alpha & mod(pvs(:,1),3)==1);
IndxPosSRejBoth1 = find(pvs(:,4)>0 & pvs(:,7) < Alpha & pvs(:,8) < Alpha & mod(pvs(:,1),3)==1);
FalseRejects1 = find(pvs(:,4)>0 & pvs(:,7) < Alpha & pvs(:,8) >= Alpha & mod(pvs(:,1),3)==1);

IndxPosSRejUnCond10 = find(pvs(:,4)>0 & pvs(:,7) < Alpha & mod(pvs(:,1),3)==2);
IndxPosSRejCond10 = find(pvs(:,4)>0 & pvs(:,8) < Alpha & mod(pvs(:,1),3)==2);
IndxPosSRejBoth10 = find(pvs(:,4)>0 & pvs(:,7) < Alpha & pvs(:,8) < Alpha & mod(pvs(:,1),3)==2);
FalseRejects10 = find(pvs(:,4)>0 & pvs(:,7) < Alpha & pvs(:,8) >= Alpha & mod(pvs(:,1),3)==2);


IndxPosSRejUnCond100 = find(pvs(:,4)>0 & pvs(:,7) < Alpha & mod(pvs(:,1),3)==0);
IndxPosSRejCond100 = find(pvs(:,4)>0 & pvs(:,8) < Alpha & mod(pvs(:,1),3)==0);
IndxPosSRejBoth100 = find(pvs(:,4)>0 & pvs(:,7) < Alpha & pvs(:,8) < Alpha & mod(pvs(:,1),3)==0);
FalseRejects100 = find(pvs(:,4)>0 & pvs(:,7) < Alpha & pvs(:,8) >= Alpha & mod(pvs(:,1),3)==0);


[length(IndxPosSRejUnCond1) length(IndxPosSRejCond1) length(IndxPosSRejBoth1) length(FalseRejects1)]
[length(IndxPosSRejUnCond10) length(IndxPosSRejCond10) length(IndxPosSRejBoth10) length(FalseRejects10)]
[length(IndxPosSRejUnCond100) length(IndxPosSRejCond100) length(IndxPosSRejBoth100) length(FalseRejects100)]

sfsFileName = '../work/Data/SFS_n_10.t_1_10_100.r_0.g_0';
sfs = dlmread(sfsFileName, ' ');
size(sfs)
[FalseRejects1 sfs(FalseRejects1,:)]
[FalseRejects10 sfs(FalseRejects10,:)]
[FalseRejects100 sfs(FalseRejects100,:)]


