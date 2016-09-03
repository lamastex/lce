%%%FileName = '../work/Tests/SFS_n_10.t_100.r_0.g_0_mcE5.pvalues.txt';

%Null Model \phi_0 = (100,0,0)
%Filename = '../work/Tests/SFS_n_10.t_100.g_0.r_0.m_100000.reps_10000.pvalues.txt'; Params = '\phi_0 = (100,0,0)';

%Alternative Models with growth g=0 and recomb rate r increases to 10 and 100
%Model 2
%Filename = '../work/Tests/SFS_n_10.t_100.g_0.r_10.m_100000.reps_10000.pvalues.txt'; Params = '\phi = (100,0,10)';
% Model 
%Filename = '../work/Tests/SFS_n_10.t_100.g_0.r_100.m_100000.reps_10000.pvalues.txt'; Params = '\phi = (100,0,100)';

%Alternative Models with recomb rate r=0 and growth g increases to 10 and 100
%
% Filename = '../work/Tests/SFS_n_10.t_100.g_10.r_0.m_100000.reps_10000.pvalues.txt'; Params = '\phi = (100,10,0)';


%Alternative Models with growth g=10 and recomb rate r increases to 10 and 100
% Model 
% Filename = '../work/Tests/SFS_n_10.t_100.g_10.r_10.m_100000.reps_10000.pvalues.txt'; Params = '\phi = (100,10,10)';
% Model 
% Filename = '../work/Tests/SFS_n_10.t_100.g_10.r_100.m_100000.reps_10000.pvalues.txt'; Params = '\phi = (100,10,100)';

%Alternative Models with growth g=100 and recomb rate r increases to 10 and 100
% Model 
% Filename = '../work/Tests/SFS_n_10.t_100.g_100.r_0.m_100000.reps_10000.pvalues.txt'; Params = '\phi = (100,100,0)';
% Filename = '../work/Tests/SFS_n_10.t_100.g_100.r_10.m_100000.reps_10000.pvalues.txt' ; Params = '\phi = (100,100,10)'
 Filename = '../work/Tests/SFS_n_10.t_100.g_100.r_100.m_100000.reps_10000.pvalues.txt'; Params = '\phi = (100,100,100)';

figure
UncondCol=7;
CondCol=8;
PlotTitle=strcat('P_{\phi_0}(|D| \leq |d_{obs}|); ',Params);
MyTajDPlots(UncondCol, CondCol, 0.0950, PlotTitle,Filename)

figure
UncondCol=11;
CondCol=12;
PlotTitle=strcat('P_{\phi_0}(D \leq d_{obs}); ',Params);
MyTajDPlots(UncondCol, CondCol, 0.2022, PlotTitle,Filename)


figure
UncondCol=13;
CondCol=14;
PlotTitle=strcat('P_{\phi_0}(D \geq d_{obs}); ',Params);
MyTajDPlots(UncondCol, CondCol, 0.0266, PlotTitle,Filename)

