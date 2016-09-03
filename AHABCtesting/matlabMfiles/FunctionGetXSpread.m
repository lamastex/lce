function [lower,upper] = FunctionGetXSpread( filenames, d1, d2 )
%Take a cell array of string filenames, open all of them, return the interval hull of the 'domains'.

%Input:
%inputfile:     a cell array of filenames
%d1:             the column index to take lower domain coordinates from
%d2:             the column index to take upper domain coordinates from
%
%Output:
%The maximum height from any of the files listed in inputfile.


n=size(filenames{1},1);
OneX1=[];
OneX2=[];
for i=1:n
    
    boxesFileName = filenames{1}{i};
    dataR = dlmread(char(boxesFileName), '\t', 0, 1); % from row 0, col 1
    
    OneX1 = [OneX1;dataR(:,d1)]; % 
    OneX2 = [OneX2;dataR(:,d2)]; % 

end
lower = min(OneX1);
upper = max(OneX2);

end

