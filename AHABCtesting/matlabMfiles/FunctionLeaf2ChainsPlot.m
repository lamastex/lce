function [ p1 p2 ] = FunctionLeaf2ChainsPlot(inputfile, axeshandle )
%Read data from inputfile and plot onto axeshandle.

%Function is set up for input file where:
%                   first column is state or line number.
%                   second column is leaf count for chain 1.
%                   fourth column is leaf count for chain 2.

%Input:
%inputfile:     full text filename of file where intput is found
%               eg 'myinput.txt'
%axeshandle:    handle of axes on which to put the plot
%
%Output:
%[p1 p2], handles of the two graphics objects, one from each chain.

dataR = dlmread(char(inputfile), '\t', 3, 0); % from row 4, col 0

states = dataR(:,1);
leaves1 = dataR(:,2);
leaves2 = dataR(:,4);

%color_edge = [.25 .25 .25];    % edge colour
axes(axeshandle);

p1 = plot(states,leaves1);
hold on
p2  = plot(states,leaves2);
hold off
legend('Chain 1', 'Chain 2');
set(get(axeshandle,'Ylabel'),'String','Leaf trace') ;
set(get(axeshandle,'Xlabel'),'String','State') ;
set(p1,'Color',[0 0 1]);
set(p2,'Color',[.7 0 .7]);
set(axeshandle,'Xlim',[min(states) max(states)]);
end

