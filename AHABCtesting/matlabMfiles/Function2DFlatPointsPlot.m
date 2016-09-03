function [ p ] = Function2DFlatPointsPlot(inputfile, axeshandle, varargin )
%Read data from inputfile and plot onto axeshandle.

%Function is set up for 2-D data.
%Plots 2-D points onto a 2-D plot

%Input:
%inputfile:     full text filename of file where intput is found
%               eg 'myinput.txt'
%axeshandle:    handle of axes on which to put the plot
%
%Up to 1 optional input:
%               colour points marker
%
%Output:
%The handle of the plot graphic.

dataR = dlmread(char(inputfile), '\t', 0, 0); % from row 0, col 0

OneX1 = dataR(:,1);

OneY1 = dataR(:,2);

Oneboxes = size(OneX1,1);

% marker colour
pcol = [0.25 0.25 1];
if  size(varargin,2) > 0
    pcol = varargin{1};
end

axes(axeshandle);

p = plot(OneX1, OneY1, 'LineStyle','none','MarkerFaceColor',pcol,'MarkerEdgeColor',pcol,'Marker','o','MarkerSize',3);


end

