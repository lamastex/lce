function [ m ] = FunctionGetMaxHeight( filenames, d )
%Take a cell array of string filenames, open all of them, return max height.

%Input:
%inputfile:     a cell array of filenames
%d:             the column index to take heights from
%
%Output:
%The maximum height from any of the files listed in inputfile.


n=size(filenames{1},1);
OneZ2=[];
for i=1:n
    
    boxesFileName = filenames{1}{i};
    dataR = dlmread(char(boxesFileName), '\t', 0, 1); % from row 0, col 1
       
    OneZ2 = [OneZ2;dataR(:,d)]; % 

end
m = max(OneZ2);

end

