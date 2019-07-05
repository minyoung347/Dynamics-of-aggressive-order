function [time_series_result] = aggregate(time_series,agg)

len=size(time_series,1);
time_series_ag=zeros(floor(len/agg),size(time_series,2));
iter=1;

for i=1:len
    if mod(i,agg)==0
        time_series_ag(iter,:)=time_series_ag(iter,:)+time_series(i,:);
        time_series_ag(iter,:)=time_series_ag(iter,:)./agg;
        iter=iter+1;
    elseif mod(i,agg)~=0 && iter<=floor(len/agg)
        if iter < floor(len/agg)
            cut = 1;
        elseif iter == floor(len/agg)
            cut = 0;
        end
        time_series_ag(iter,:)=time_series_ag(iter,:)+time_series(i,:);
    end
end

time_series_result=time_series_ag(1:end-cut,:);


