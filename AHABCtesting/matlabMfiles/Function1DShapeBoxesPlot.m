function [ patchhandle ] = Function1DShapeBoxesPlot(inputfile, axeshandle, varargin )
%Read data from inputfile and plot onto axeshandle.

%Function is set up for 1-D data from a shape, 
% ie vertex input, not interval ends.
%Plots 1-D box and real range as a 2-d box

%Input:
%inputfile:     full text filename of file where intput is found
%               eg 'myinput.txt'
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

OneZ2 = dataR(:,1); % 
OneZ1 = zeros(size(OneZ2,1));

vert1 = dataR(:,2); 

vert2 = dataR(:,3); 

Oneboxes = size(OneZ2,1);

OneVert=[vert1(1) OneZ1(1);... % 1
    vert2(1) OneZ1(1);...   % 2
    vert2(1) OneZ2(1);...   % 3
    vert1(1) OneZ2(1)];     % 4

OneFaceBase = [1 2 3 4];               
   


OneFace = OneFaceBase;

fcol = [.7 0.6 1];
if  size(varargin,2) > 0
    fcol = varargin{1};
end

tcolbase = fcol;
tcol = tcolbase;

for i=2:Oneboxes
    new = [vert1(i) OneZ1(i);... % 1
    vert2(i) OneZ1(i);...   % 2
    vert2(i) OneZ2(i);...   % 3
    vert1(i) OneZ2(i)];     % 4

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
xlim([min([vert1; vert2]) max([vert1; vert2])]);
ylim([0 max(OneZ2)]);


end

