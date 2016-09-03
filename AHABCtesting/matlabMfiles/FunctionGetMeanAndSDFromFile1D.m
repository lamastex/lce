function [ mn sd ] = FunctionGetMeanAndSDFromFile1D(inputfile )
%Read data from inputfile and return mean and sd

%Function is set up for 1-D data 

%Input:
%inputfile:     full text filename of file where intput is found
%               eg 'myinput.txt'
%
%Output:
%The mean and sd.

dataR = dlmread(char(inputfile), '\t', 0, 0); % from row 0, col 0

stats = dataR(:,1);

mn = mean(stats);
sd = std(stats);


end

