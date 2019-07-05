function [B] = burstiness(time_series,alt)
% burstness of inter-event time

r = std(time_series) / mean(time_series);
n = length(time_series);
y_tildar = min(time_series) / sum(time_series);

if alt == 0
    B = (r-1) / (r+1);
elseif alt == 1
    B = (sqrt(n+1) * r - sqrt(n-1)) / ((sqrt(n+1) - 2) * r + sqrt(n-1));
elseif alt == 2
    B = (n-2)*(sqrt(n+1) * r - (sqrt(n-1))*(1-n*y_tildar)) / ((n*sqrt(n+1) - 2*(n-1)) * r + sqrt(n-1)*(n-2*sqrt(n+1))*(1-n*y_tildar));
else
    disp('Error: alt is 0 or 1 or 2')
end

% B = (std(time_series)-mean(time_series))/(std(time_series)+mean(time_series));