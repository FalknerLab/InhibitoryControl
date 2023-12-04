function [val] = sem(data_set,dim)
%function to determine the standard error of the mean for a data set
%created 08-17-16, modified 04-09-18

if nargin == 1
    error('need to specify direction')
elseif nargin == 2
    if dim == 1 %calc sem for columns
        val = std(data_set)/sqrt(size(data_set,1));
    elseif dim == 2 %calc sem for rows
        val = std(data_set,0,2)/sqrt(size(data_set,2));
    end
end
end