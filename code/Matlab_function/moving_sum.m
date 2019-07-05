function [x_result] = moving_sum(x,win_size,mv_size,start_end)
% Calculate Moving Average
% column of x is time axis, row of x is item axis
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

x_result = zeros(ceil((size(x,1)-win_size+1)/mv_size),size(x,2));

if start_end == 1    
    for i=1:ceil((size(x,1)-win_size+1)/mv_size)
        x_result(i,:) = sum(x(1+(i-1)*mv_size:1+(i-1)*mv_size+win_size-1,:),1);
    end
elseif start_end == 2
    for i=1:ceil((size(x,1)-win_size+1)/mv_size)
        x_result(end-i+1,:) = sum(x(end-(i-1)*mv_size-win_size+1:end-(i-1)*mv_size,:),1);
    end
end



% x_result = zeros(floor(size(x,1)/win_size),size(x,2));
% if se_ind ==1
%     for i=1:floor(size(x,1)/win_size)
%         x_result(i,:) = mean(x(1+(i-1)*win_size:i*win_size,:),1);
%     end
% elseif se_ind ==2
%     for i=1:floor(size(x,1)/win_size)
%         x_result(end-i+1,:) = mean(x(end-i*win_size+1:end-(i-1)*win_size,:),1);
%     end
% end


% x_result = zeros(size(x,1)-win_size+1,size(x,2));
% for i=1:size(x,1)-win_size
%     x_result(i,:) = mean(x(i:i+win_size,:),1);
% end