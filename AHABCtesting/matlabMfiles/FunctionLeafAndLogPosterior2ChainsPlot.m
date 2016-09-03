function [ AX, H11, H12 H21, H22] = FunctionLeafAndLogPosterior2ChainsPlot(inputfile, axeshandle )
%Read data from inputfile and plot onto axeshandle.

%Function is set up for input file where:
%                   first column is state or line number.
%                   second column is leaf count for chain 1
%                   third column is log posterior for chain 1
%                   fourth column is leaf count for chain 2
%                   fifth column is log posterior for chain 2

%Input:
%inputfile:     full text filename of file where intput is found
%               eg 'myinput.txt'
%axeshandle:    handle of axes on which to put the plot
%
%Output:
%[AX,H1,H2], handles of the two axes created in axeshandles and the handles
%           of the graphics objects from each plot in H1 and H2.
%           AX(1) is the left axes and AX(2) is the right axes.

dataR = dlmread(char(inputfile), '\t', 3, 0); % from row 4, col 0

states = dataR(:,1);
leaves1 = dataR(:,2);
logposteriors1 = dataR(:,3);
leaves2 = dataR(:,4);
logposteriors2 = dataR(:,5);

%color_edge = [.25 .25 .25];    % edge colour
axes(axeshandle);

[ AX, H11, H12 ] = plotyy(states,leaves1, states, logposteriors1);
hold on
[ AX, H21, H22 ]  = plotyy(states,leaves2, states, logposteriors2);
hold off
set(get(AX(1),'Ylabel'),'String','Leaf trace') ;
set(get(AX(2),'Ylabel'),'String','Log-posterior') ;
set(get(AX(1),'Xlabel'),'String','State') ;
set(H21,'Color',[.7 0 .7])
set(H22,'Color',[0 .7 .7])
end

