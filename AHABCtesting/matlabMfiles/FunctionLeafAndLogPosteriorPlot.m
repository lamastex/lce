function [ AX, H1, H2 ] = FunctionLeafAndLogPosteriorPlot(inputfile, colindex, axeshandle )
%Read data from inputfile and plot onto axeshandle.

%Function is set up for input file where:
%                   first column is state or line number.
%                   colindex column is leaf count
%                   colindex+1 column is log posterior

%Input:
%inputfile:     full text filename of file where intput is found
%               eg 'myinput.txt'
%colindex:      colindex where we expect to find the leaf counts
%axeshandle:    handle of axes on which to put the plot
%
%Output:
%[AX,H1,H2], handles of the two axes created in axeshandles and the handles
%           of the graphics objects from each plot in H1 and H2.
%           AX(1) is the left axes and AX(2) is the right axes.

dataR = dlmread(char(inputfile), '\t', 3, 0); % from row 3, col 0

states = dataR(:,1);
leaves = dataR(:,colindex);
logposteriors = dataR(:,colindex+1);

%color_edge = [.25 .25 .25];    % edge colour
axes(axeshandle);

[ AX, H1, H2 ] = plotyy(states,leaves, states, logposteriors);
set(get(AX(1),'Ylabel'),'String','Leaf trace') ;
set(get(AX(2),'Ylabel'),'String','Log-posterior') ;
set(get(AX(1),'Xlabel'),'String','State') ;

end

