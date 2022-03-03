function c = numlist2cellstr(numlist)
% numlist2cellstr is to make a cell array (especially for 
%   a string list for popupmenu) from a number array, numlist
% try numlist2cellstr(1:5);
%
%    cellstr(num2str(11))
c = cell(numel(numlist), 1);
for i=1:numel(numlist)
    c{i,1} = num2str(i);
end