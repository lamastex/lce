function [ f v1 ] = ReadFilenamesAndVal( inputfile )
%Read string filenames from a file and return as a cell array of strings.

%Input:
%inputfile:     full text filename of file where filenames are
%               eg 'myinput.txt'
%
%Output:
%The cell array of strings read from the inputfile
%A matrix of values read from the input file.

fid = fopen(inputfile);
C = textscan(fid, '%s %f','EndOfLine','\n');
fclose(fid);
f = C(:,1);
v1 = cell2mat(C(:,2));

end

