function [result, ind_surv, ind_del] = delniz(data,flag,surv_on)
% reshape data that excluding 
% [result, ind_surv, ind_del] = delniz(data,flag,surv_on)
% Nan if flag == 1
% Inf if flag == 2
% Zero if flag == 3
% data is n x 1 array
% flag is an array whose length is 1 to 3
% ind_surv is an array that survive from the filter
% surv_on is on(1)/off(2) for ind_surv
% ind_del is an array deleted by the filter
% 2019.03.22, Min-Young Lee
% minyoung@postech.ac.kr

if size(data,2) > size(data,1)
    data = data';
end

ind_surv = (1:length(data))';
ind_del = [];
if sum(find(flag == 1)) > 0
    ind_del = [ind_del; find(isnan(data))];
end
if sum(find(flag == 2)) > 0 
    ind_del = [ind_del; find(isinf(data))];
end
if sum(find(flag == 3)) > 0
    ind_del = [ind_del; find(data == 0)];
end
ind_del = sort(ind_del);

if surv_on == 1
    for i = 1:length(ind_del)
        ind_surv(ind_surv == ind_del(i)) = [];
    end
end

for i = 1:length(flag)
    if flag(i) == 1
        data = data(~isnan(data));
    elseif flag(i) == 2
        data = data(~isinf(data));
    elseif flag(i) == 3
        data = data(~isinf(1./data));
    end
end

result = data;

