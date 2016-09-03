function [ patchhandle ] = Function1DShapeFromPavingBoxesPlot(inputfile, scale, shift, axeshandle, varargin )
%Read data from inputfile and plot onto axeshandle.

%Function is set up for 1-D data from a paving, to be transformed into a shape, 
% ie interval end input (transformed into vertex data here).
%Plots 1-D shape and real range as a 2-d shape

%Input:
%inputfile:     full text filename of file where intput is found
%               eg 'myinput.txt'
%scale:         a 1-d vector of scales (to scale up the shape)
%shift:         a 1-d vector of shifts (to move the shape)
%axeshandle:    handle of axes on which to put the plot
%
%Up to 4 optional input:
%               colour for top and bottom faces of the boxes plotted
%               colour for edges of the faces of the boxes plotted
%               face alpha (between 0 for clear, 1 opaque, inclusive)
%               edge alpha (between 0 for clear, 1 opaque, inclusive)
%
%Output:
%The handle of the patch graphic.

dataR = dlmread(char(inputfile), '\t', 0, 1); % from row 0, col 1

OneZ2 = dataR(:,2)/(scale(1)*scale(2));
OneZ1 = zeros(size(OneZ2,1));

vert11 = dataR(:,3)*scale(1) + shift(1);

vert21 = dataR(:,4)*scale(1) + shift(1);

Oneboxes = size(OneZ2,1);

OneVert=[vert11(1) OneZ1(1);... % 1
    vert21(1) OneZ2(1);...   % 2
    vert11(1) OneZ1(1);...   % 4
    vert21(1) OneZ2(1)];     % 7

OneFaceBase = [1 2 4 3];              
   


OneFace = OneFaceBase;

fcol = [.7 0.6 1];
if  size(varargin,2) > 0
    fcol = varargin{1};
end

tcolbase = fcol;
tcol = tcolbase;

for i=2:Oneboxes
    new = [vert11(i) OneZ1(i);... % 1
    vert21(i) OneZ2(i);...   % 2
    vert11(i) OneZ1(i);...   % 4
    vert21(i) OneZ2(i)];     % 7

    OneVert = [OneVert; new];

    f = OneFaceBase+(4*(i-1));
    OneFace = [OneFace; f];

    t=tcolbase;
    tcol = [tcol; t];

end

ecol = [.25 .25 .25];    % edge colour
if  size(varargin,2) > 1
    ecol = varargin{2};
end

falpha = 0.6;    
if  size(varargin,2) > 2
    falpha = varargin{3};
end

ealpha = 1;   
if  size(varargin,2) > 3
    ealpha = varargin{4};
end


axes(axeshandle);

patchhandle = patch('Faces',OneFace,'Vertices',OneVert,'FaceVertexCData',tcol,...
      'FaceColor','flat','FaceAlpha', falpha,'EdgeColor', ecol,'EdgeAlpha',ealpha);
%xlim([min([vert11; vert21]) max([vert11; vert21])]);
%ylim([0 max(OneZ2)]);

end

