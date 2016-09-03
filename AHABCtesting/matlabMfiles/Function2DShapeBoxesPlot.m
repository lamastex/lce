function [ patchhandle ] = Function2DShapeBoxesPlot(inputfile, axeshandle, varargin )
%Read data from inputfile and plot onto axeshandle.

%Function is set up for 2-D data from a shape, 
% ie vertex input, not interval ends.
%Plots 2-D box and real range as a 3-d box

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

OneZ2 = dataR(:,1);
OneZ1 = zeros(size(OneZ2,1));

vert11 = dataR(:,2);
vert12 = dataR(:,3);

vert21 = dataR(:,4);
vert22 = dataR(:,5);

vert31 = dataR(:,6);
vert32 = dataR(:,7);

vert41 = dataR(:,8);
vert42 = dataR(:,9);

Oneboxes = size(OneZ2,1);

OneVert=[vert11(1) vert12(1) OneZ1(1);... % 1
    vert21(1) vert22(1) OneZ1(1);...   % 2
    vert31(1) vert32(1) OneZ1(1);...   % 4
    vert41(1) vert42(1) OneZ1(1);...   % 3
    vert11(1) vert12(1) OneZ2(1);... % 5
    vert21(1) vert22(1) OneZ2(1);...   % 6
    vert31(1) vert32(1) OneZ2(1);...   % 8
    vert41(1) vert42(1) OneZ2(1)];     % 7

OneFaceBase = [1 2 4 3;...     % bottom
    5 6 8 7;...             % top
    1 2 6 5;...             % front
    4 3 7 8;...             % back
    2 4 8 6;...             % rhs    
    3 1 5 7];               % lhs
   


OneFace = OneFaceBase;

fcol = [.7 0.6 1];
if  size(varargin,2) > 0
    fcol = varargin{1};
end

tcolbase = [fcol;...
    fcol;...
    1 1 1;...
    1 1 1;...
    1 1 1;...
    1 1 1];
tcol = tcolbase;

for i=2:Oneboxes
    new = [vert11(i) vert12(i) OneZ1(i);... % 1
    vert21(i) vert22(i) OneZ1(i);...   % 2
    vert31(i) vert32(i) OneZ1(i);...   % 4
    vert41(i) vert42(i) OneZ1(i);...   % 3
    vert11(i) vert12(i) OneZ2(i);... % 5
    vert21(i) vert22(i) OneZ2(i);...   % 6
    vert31(i) vert32(i) OneZ2(i);...   % 8
    vert41(i) vert42(i) OneZ2(i)];     % 7

    OneVert = [OneVert; new];

    f = OneFaceBase+(8*(i-1));
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
xlim([min([vert11; vert21; vert31; vert41]) max([vert11; vert21; vert31; vert41])]);
ylim([min([vert12; vert22; vert32; vert42]) max([vert12; vert22; vert32; vert42])]);
zlim([0 max(OneZ2)]);

end

