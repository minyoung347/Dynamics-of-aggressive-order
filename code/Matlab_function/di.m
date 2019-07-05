function [result_series] = di(series)

result_series=series(2:end,:)-series(1:end-1,:);