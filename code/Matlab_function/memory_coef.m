function [M] = memory_coef(time_series)
% Memory coefficient of inter-event times

len = length(time_series);
x1 = time_series(1:end-1);
x2 = time_series(2:end);

M = 1/(len-1)/(std(x1)*std(x2)) * sum((x1-mean(x1)).*(x2-mean(x2)));