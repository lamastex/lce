function [ stats ] = FunctionReadMeansAndSds(inputfile)
%Read sds and means from inputfile.

%Function is set up for file with sds in first column, means in second
%Input:
%inputfile:     full text filename of file where intput is found
%               eg 'myinput.txt'
%
%Output:
%A (2 x n) matrix of stats, sds first row,  means in second row, where n
% is the number of means = number of sds, ie the number  rows in the input.

dataR = dlmread(char(inputfile), '\t', 0, 0); % from row 0, col 0

sds = dataR(:,1);
means = dataR(:,2);

stats = [sds';means'];

end

