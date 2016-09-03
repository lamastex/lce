% script to plot one 2d function using real ranges
% 
% 

%put the name of the input file here
boxesFileName = 'simdata4DGR_Likelihood_s_10_ns_25000_nr_1000_th_0.02_g_50.0.txt';

%change the figure handle if necessary
figure;

clear functions
f1 = @Function2DBoxesPlot;

h1 = gca;
cla(h1);

p = f1(boxesFileName, h1);
set(get(h1,'XLabel'),'String',texlabel('theta_1'));
set(get(h1,'YLabel'),'String',texlabel('theta_2'));
%set(h1,'View',[61 18]); % good for untransformed data
set(h1,'View',[21 32]); 
da = daspect;
if da(1) > da(2)
    da(2) = da(1);
elseif da(2) > da(1)
    da(1) = da(2); 
end

daspect(h1, da);
set(get(h1,'Title'),'String',boxesFileName);
