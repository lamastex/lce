function [ f s1 ] = ReadFilenamesAndSize( inputfile )
%Read string filenames from a file and return as a cell array of strings.

%Input:
%inputfile:     full text filename of file where filenames are
%               eg 'myinput.txt'
%
%Output:
%The cell array of strings read from the inputfile.

fid = fopen(inputfile);
C = textscan(fid, '%s %u','EndOfLine','\n');
fclose(fid);
f = C(:,1);
s1 = cell2mat(C(:,2));

end

