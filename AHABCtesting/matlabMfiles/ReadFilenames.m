function [ C ] = ReadFilenames( inputfile )
%Read string filenames from a file and return as a cell array of strings.

%Input:
%inputfile:     full text filename of file where filenames are
%               eg 'myinput.txt'
%
%Output:
%The cell array of strings read from the inputfile.

fid = fopen(inputfile);
C = textscan(fid, '%s','EndOfLine','\n');
fclose(fid);

end

