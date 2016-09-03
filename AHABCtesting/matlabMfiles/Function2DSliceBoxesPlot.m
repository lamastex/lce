function [ patchhandle ] = Function2DSliceBoxesPlot(inputfile, axeshandle )
%Read data from inputfile and plot onto axeshandle.

%Function is set up for 2-D data with real range.
%Plots 2-D box and real range as a 3-d box

%Input:
%inputfile:     full text filename of file where intput is found
%               eg 'myinput.txt'
%axeshandle:    handle of axes on which to put the plot
%
%Output:
%The handle of the patch graphic.

dataR = dlmread(char(inputfile), '\t', 0, 1); % from row 0, col 1

OneZ2 = dataR(:,2);
OneZ1 = zeros(size(OneZ2,1));

OneX1 = dataR(:,3);
OneX2 = dataR(:,4);

OneY1 = dataR(:,5);
OneY2 = dataR(:,6);

Oneboxes = size(OneX1,1);

OneVert=[OneX1(1) OneY1(1) OneZ1(1);... % 1
    OneX2(1) OneY1(1) OneZ1(1);...   % 2
    OneX1(1) OneY2(1) OneZ1(1);...   % 3 4
    OneX2(1) OneY2(1) OneZ1(1);...   % 4 3
    OneX1(1) OneY1(1) OneZ2(1);...   % 5
    OneX2(1) OneY1(1) OneZ2(1);...   % 6
    OneX1(1) OneY2(1) OneZ2(1);...   % 7 8
    OneX2(1) OneY2(1) OneZ2(1)];     % 8 7

OneFaceBase = [1 2 4 3;...     % bottom
    5 6 8 7;...             % top
    1 2 6 5;...             % front
    4 3 7 8;...             % back
    2 4 8 6;...             % rhs    
    3 1 5 7];               % lhs
   


OneFace = OneFaceBase;

tcolbase = [1 0 0;...
    1 0 0;...
    1 0 0;...
    1 0 0;...
    1 0 0;...
    1 0 0];
tcol = tcolbase;

for i=2:Oneboxes
    new = [OneX1(i) OneY1(i) OneZ1(i);... % 1
    OneX2(i) OneY1(i) OneZ1(i);...   % 2
    OneX1(i) OneY2(i) OneZ1(i);...   % 3
    OneX2(i) OneY2(i) OneZ1(i);...   % 4
    OneX1(i) OneY1(i) OneZ2(i);...   % 5
    OneX2(i) OneY1(i) OneZ2(i);...   % 6
    OneX1(i) OneY2(i) OneZ2(i);...   % 7
    OneX2(i) OneY2(i) OneZ2(i)];     % 8

    OneVert = [OneVert; new];

    f = OneFaceBase+(8*(i-1));
    OneFace = [OneFace; f];

    t=tcolbase;
    tcol = [tcol; t];

end

color_edge = [.25 .25 .25];    % edge colour
axes(axeshandle);

patchhandle = patch('Faces',OneFace,'Vertices',OneVert,'FaceVertexCData',tcol,...
      'FaceColor','flat','FaceAlpha', 1,'EdgeColor', color_edge,'EdgeAlpha',1);


end

