function [ patchhandle ] = Function2DFlatBoxesPlot(inputfile, axeshandle, varargin )
%Read data from inputfile and plot onto axeshandle.

%Function is set up for 2-D data.
%Plots 2-D box (no range)

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

OneX1 = dataR(:,3);
OneX2 = dataR(:,4);

OneY1 = dataR(:,5);
OneY2 = dataR(:,6);

Oneboxes = size(OneX1,1);

OneVert=[OneX1(1) OneY1(1);... % 1
    OneX2(1) OneY1(1);...   % 2
    OneX1(1) OneY2(1);...   % 3 4
    OneX2(1) OneY2(1)];...   % 4 3
    

OneFaceBase = [1 2 4 3];
    
   


OneFace = OneFaceBase;

fcol = [1 1 1];
if  size(varargin,2) > 0
    fcol = varargin{1};
end

tcolbase = fcol;
tcol = tcolbase;


for i=2:Oneboxes
    new = [OneX1(i) OneY1(i);... % 1
    OneX2(i) OneY1(i);...   % 2
    OneX1(i) OneY2(i);...   % 3
    OneX2(i) OneY2(i)];     % 8

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
xlim([min(OneX1) max(OneX2)]);
ylim([min(OneY1) max(OneY2)]);


end

