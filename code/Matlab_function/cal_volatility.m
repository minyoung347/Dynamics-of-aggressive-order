function volatility = cal_volatility(x,win_size,mv_size,start_end)
% Calculate volatility of price return
% x is [x by 1] price return vector
% win_size is the size of window
% mv_size is the moving distance
% se_ind refers conserving start point or end point
% win_size < size(x,1)
% mv_size <= win_size
% size(x_result,1) = ceil( (size(x,1) - win_size + 1) / mv_size )
% if mv_size == win_size, size(x_result,1) = floor(size(x,1) / mv_size)
% if mv_size == 1, size(x_result,1) = size(x,1) - win_size + 1
% if start_end ==1, conserve start points and throw away end points
% if start_end ==2, conserve end points and throw away start points


volatility = zeros(ceil((size(x,1)-win_size+1)/mv_size),1);

if start_end == 1
    for i=1:ceil((size(x,1)-win_size+1)/mv_size)
        volatility(i,:) = std(x(1+(i-1)*mv_size:1+(i-1)*mv_size+win_size-1,:),1);
    end
elseif start_end == 2
    for i=1:ceil((size(x,1)-win_size+1)/mv_size)
        volatility(end-i+1,:) = std(x(end-(i-1)*mv_size-win_size+1:end-(i-1)*mv_size,:),1);
    end
end
