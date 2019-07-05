function [slope] = dist_slope(x,y)

slope = diff(y)./diff(x);