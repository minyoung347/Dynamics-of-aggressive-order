function [spread_plus,spread_minus,bid_vol_minus,bid_vol_plus,ask_vol_minus,ask_vol_plus,bid_vol_price_minus,bid_vol_price_plus,ask_vol_price_minus,ask_vol_price_plus,bid_gap_minus,bid_gap_plus,ask_gap_minus,ask_gap_plus,bid_price_minus,bid_price_plus,ask_price_minus,ask_price_plus] = calculating_event_ver6_1905(t,trade,spread_a,bid_vol,ask_vol,bid_vol_price,ask_vol_price,bid_gap,ask_gap,bid_price,ask_price,tick_bin_valid,ind_minus,ind_plus,NI)
% vi = [10, 100, 500, 1000, 3000, 5000, 8000, 10000, 20000, 30000];

for ti=1:length(tick_bin_valid)
    for i=1:NI
        ind_m=eval(sprintf('ind_minus{%d,1}{%d,1}',i,ti)); % trading in the bid side (sell initiated market order)
        ind_p=eval(sprintf('ind_plus{%d,1}{%d,1}',i,ti)); % trading in the ask side (buy initiated market order)
        if isempty(ind_m)==0
            for ii=1:length(ind_m)
                ind_mm{ti,1}{i,1}(ii,1)=find(trade(:,4)==ind_m(ii,1));
            end
        else
            ind_mm{ti,1}{i,1}=[];
        end
        if isempty(ind_p)==0
            for ii=1:length(ind_p)
                ind_pp{ti,1}{i,1}(ii,1)=find(trade(:,4)==ind_p(ii,1));
            end
        else
            ind_pp{ti,1}{i,1}=[];
        end
        clear ind_m ind_p
    end
end

% cut the length to calculate the time delayed mean
len=size(trade,1);
for ti=1:length(tick_bin_valid)
    for i=1:NI
        ind_mm{ti,1}{i,1}=ind_mm{ti,1}{i,1}(ind_mm{ti,1}{i,1}>t+1 & ind_mm{ti,1}{i,1}<len-t-1);
        ind_pp{ti,1}{i,1}=ind_pp{ti,1}{i,1}(ind_pp{ti,1}{i,1}>t+1 & ind_pp{ti,1}{i,1}<len-t-1);
    end
end
% cut the negative spread
for ti=1:length(tick_bin_valid)
    for i=1:NI
        for ii=1:2*t+1
            lag=-t-1+ii;
            if isempty(ind_pp{ti,1}{i,1})==0
                ind_ppp{ti,1}{i,ii}=ind_pp{ti,1}{i,1}(find(spread_a(ind_pp{ti,1}{i,1}+lag)>0),1);
            else
                ind_ppp{ti,1}{i,ii}=[];
            end
            if isempty(ind_mm{ti,1}{i,1})==0
                ind_mmm{ti,1}{i,ii}=ind_mm{ti,1}{i,1}(find(spread_a(ind_mm{ti,1}{i,1}+lag)>0),1);
            else
                ind_mmm{ti,1}{i,ii}=[];
            end
        end
    end              
end
% cut the negative spread and delete the overlapping data
for ti=1:length(tick_bin_valid)
    for i=1:NI
        if isempty(ind_mmm{ti,1}{i,1})==0
            temp1=ind_mmm{ti,1}{i,1};
            for ii=1:2*t
                temp1=intersect(temp1,ind_mmm{ti,1}{i,ii+1});
            end
        else
            temp1=[];
        end
        if isempty(ind_ppp{ti,1}{i,1})==0
            temp2=ind_ppp{ti,1}{i,1};
            for ii=1:2*t
                temp2=intersect(temp2,ind_ppp{ti,1}{i,ii+1});
            end
        else
            temp2=[];
        end
        ind_mmmm{ti,1}{i,1}=temp1;
        ind_pppp{ti,1}{i,1}=temp2;
    end                   
end

for ti=1:length(tick_bin_valid)
    for i=1:NI
        for ii=1:2*t+1
            lag=-t-1+ii;
            for iii=1:10
                bid_vol_m{ti,1}{i,ii}(:,iii)=bid_vol(ind_mmmm{ti,1}{i,1}+lag,iii);
                bid_vol_p{ti,1}{i,ii}(:,iii)=bid_vol(ind_pppp{ti,1}{i,1}+lag,iii);
                ask_vol_m{ti,1}{i,ii}(:,iii)=ask_vol(ind_mmmm{ti,1}{i,1}+lag,iii);
                ask_vol_p{ti,1}{i,ii}(:,iii)=ask_vol(ind_pppp{ti,1}{i,1}+lag,iii);
                bid_vol_price_m{ti,1}{i,ii}(:,iii)=bid_vol_price(ind_mmmm{ti,1}{i,1}+lag,iii);
                bid_vol_price_p{ti,1}{i,ii}(:,iii)=bid_vol_price(ind_pppp{ti,1}{i,1}+lag,iii);
                ask_vol_price_m{ti,1}{i,ii}(:,iii)=ask_vol_price(ind_mmmm{ti,1}{i,1}+lag,iii);
                ask_vol_price_p{ti,1}{i,ii}(:,iii)=ask_vol_price(ind_pppp{ti,1}{i,1}+lag,iii);
                bid_gap_m{ti,1}{i,ii}(:,iii)=bid_gap(ind_mmmm{ti,1}{i,1}+lag,iii);
                bid_gap_p{ti,1}{i,ii}(:,iii)=bid_gap(ind_pppp{ti,1}{i,1}+lag,iii);
                ask_gap_m{ti,1}{i,ii}(:,iii)=ask_gap(ind_mmmm{ti,1}{i,1}+lag,iii);
                ask_gap_p{ti,1}{i,ii}(:,iii)=ask_gap(ind_pppp{ti,1}{i,1}+lag,iii);
            end
            spread_m{ti,1}{i,ii}=spread_a(ind_mmmm{ti,1}{i,1}+lag,1) / tick_bin_valid(ti);
            spread_p{ti,1}{i,ii}=spread_a(ind_pppp{ti,1}{i,1}+lag,1) / tick_bin_valid(ti);
            
            bid_price_m{ti,1}{i,ii}=bid_price(ind_mmmm{ti,1}{i,1}+lag,1);
            bid_price_p{ti,1}{i,ii}=bid_price(ind_pppp{ti,1}{i,1}+lag,1);
            ask_price_m{ti,1}{i,ii}=ask_price(ind_mmmm{ti,1}{i,1}+lag,1);
            ask_price_p{ti,1}{i,ii}=ask_price(ind_pppp{ti,1}{i,1}+lag,1);
        end
    end
end

for ti=1:length(tick_bin_valid)
    for i=1:NI % aggressiveness
        for ii=1:2*t+1 % time lag
            if isempty(spread_m{ti,1}{i,ii})==0
%                 flow_buy{ti,1}{i,1}(ii,1)=mean(flow_bb{ti,1}{i,ii});
%                 flow_buy{ti,1}{i,1}(ii,2)=size(flow_bb{ti,1}{i,ii},1);
                spread_minus{ti,1}{i,1}(ii,1)=mean(spread_m{ti,1}{i,ii});
                spread_minus{ti,1}{i,1}(ii,2)=size(spread_m{ti,1}{i,ii},1);
                bid_price_minus{ti,1}{i,1}(ii,1)=mean(bid_price_m{ti,1}{i,ii});
                bid_price_minus{ti,1}{i,1}(ii,2)=size(bid_price_m{ti,1}{i,ii},1);
                ask_price_minus{ti,1}{i,1}(ii,1)=mean(ask_price_m{ti,1}{i,ii});
                ask_price_minus{ti,1}{i,1}(ii,2)=size(ask_price_m{ti,1}{i,ii},1);
            else
%                 flow_buy{ti,1}{i,1}(ii,1)=0;
%                 flow_buy{ti,1}{i,1}(ii,2)=0;
                spread_minus{ti,1}{i,1}(ii,1)=0;
                spread_minus{ti,1}{i,1}(ii,2)=0;
                bid_price_minus{ti,1}{i,1}(ii,1)=0;
                bid_price_minus{ti,1}{i,1}(ii,2)=0;
                ask_price_minus{ti,1}{i,1}(ii,1)=0;
                ask_price_minus{ti,1}{i,1}(ii,2)=0;
            end
            if isempty(spread_p{ti,1}{i,ii})==0
%                 flow_sell{ti,1}{i,1}(ii,1)=mean(flow_aa{ti,1}{i,ii});
%                 flow_sell{ti,1}{i,1}(ii,2)=size(flow_aa{ti,1}{i,ii},1);
                spread_plus{ti,1}{i,1}(ii,1)=mean(spread_p{ti,1}{i,ii});
                spread_plus{ti,1}{i,1}(ii,2)=size(spread_p{ti,1}{i,ii},1);
                bid_price_plus{ti,1}{i,1}(ii,1)=mean(bid_price_p{ti,1}{i,ii});
                bid_price_plus{ti,1}{i,1}(ii,2)=size(bid_price_p{ti,1}{i,ii},1);
                ask_price_plus{ti,1}{i,1}(ii,1)=mean(ask_price_p{ti,1}{i,ii});
                ask_price_plus{ti,1}{i,1}(ii,2)=size(ask_price_p{ti,1}{i,ii},1);
            else
%                 flow_sell{ti,1}{i,1}(ii,1)=0;
%                 flow_sell{ti,1}{i,1}(ii,2)=0;
                spread_plus{ti,1}{i,1}(ii,1)=0;
                spread_plus{ti,1}{i,1}(ii,2)=0;
                bid_price_plus{ti,1}{i,1}(ii,1)=0;
                bid_price_plus{ti,1}{i,1}(ii,2)=0;
                ask_price_plus{ti,1}{i,1}(ii,1)=0;
                ask_price_plus{ti,1}{i,1}(ii,2)=0;
            end
        end
    end
end

for ti=1:length(tick_bin_valid)
    for i=1:NI
        for ii=1:2*t+1
            for iii=1:10 % order book depth
                % bid_vol
                if isempty(bid_vol_m{ti,1}{i,ii})==0
                    bid_vol_minus{ti,1}{i,iii}(ii,1)=mean(bid_vol_m{ti,1}{i,ii}(:,iii));
                    bid_vol_minus{ti,1}{i,iii}(ii,2)=size(bid_vol_m{ti,1}{i,ii}(:,iii),1);
                else
                    bid_vol_minus{ti,1}{i,iii}(ii,1)=0;
                    bid_vol_minus{ti,1}{i,iii}(ii,2)=0;
                end
                if isempty(bid_vol_p{ti,1}{i,ii})==0
                    bid_vol_plus{ti,1}{i,iii}(ii,1)=mean(bid_vol_p{ti,1}{i,ii}(:,iii));
                    bid_vol_plus{ti,1}{i,iii}(ii,2)=size(bid_vol_p{ti,1}{i,ii}(:,iii),1);
                else
                    bid_vol_plus{ti,1}{i,iii}(ii,1)=0;
                    bid_vol_plus{ti,1}{i,iii}(ii,2)=0;
                end

                % ask_vol
                if isempty(ask_vol_m{ti,1}{i,ii})==0
                    ask_vol_minus{ti,1}{i,iii}(ii,1)=mean(ask_vol_m{ti,1}{i,ii}(:,iii));
                    ask_vol_minus{ti,1}{i,iii}(ii,2)=size(ask_vol_m{ti,1}{i,ii}(:,iii),1);
                else
                    ask_vol_minus{ti,1}{i,iii}(ii,1)=0;
                    ask_vol_minus{ti,1}{i,iii}(ii,2)=0;
                end
                if isempty(ask_vol_p{ti,1}{i,ii})==0
                    ask_vol_plus{ti,1}{i,iii}(ii,1)=mean(ask_vol_p{ti,1}{i,ii}(:,iii));
                    ask_vol_plus{ti,1}{i,iii}(ii,2)=size(ask_vol_p{ti,1}{i,ii}(:,iii),1);
                else
                    ask_vol_plus{ti,1}{i,iii}(ii,1)=0;
                    ask_vol_plus{ti,1}{i,iii}(ii,2)=0;
                end

                 % bid_vol_price
                if isempty(bid_vol_price_m{ti,1}{i,ii})==0
                    bid_vol_price_minus{ti,1}{i,iii}(ii,1)=mean(bid_vol_price_m{ti,1}{i,ii}(:,iii));
                    bid_vol_price_minus{ti,1}{i,iii}(ii,2)=size(bid_vol_price_m{ti,1}{i,ii}(:,iii),1);
                else
                    bid_vol_price_minus{ti,1}{i,iii}(ii,1)=0;
                    bid_vol_price_minus{ti,1}{i,iii}(ii,2)=0;
                end
                if isempty(bid_vol_price_p{ti,1}{i,ii})==0
                    bid_vol_price_plus{ti,1}{i,iii}(ii,1)=mean(bid_vol_price_p{ti,1}{i,ii}(:,iii));
                    bid_vol_price_plus{ti,1}{i,iii}(ii,2)=size(bid_vol_price_p{ti,1}{i,ii}(:,iii),1);
                else
                    bid_vol_price_plus{ti,1}{i,iii}(ii,1)=0;
                    bid_vol_price_plus{ti,1}{i,iii}(ii,2)=0;
                end

                % ask_vol_price
                if isempty(ask_vol_price_m{ti,1}{i,ii})==0
                    ask_vol_price_minus{ti,1}{i,iii}(ii,1)=mean(ask_vol_price_m{ti,1}{i,ii}(:,iii));
                    ask_vol_price_minus{ti,1}{i,iii}(ii,2)=size(ask_vol_price_m{ti,1}{i,ii}(:,iii),1);
                else
                    ask_vol_price_minus{ti,1}{i,iii}(ii,1)=0;
                    ask_vol_price_minus{ti,1}{i,iii}(ii,2)=0;
                end
                if isempty(ask_vol_price_p{ti,1}{i,ii})==0
                    ask_vol_price_plus{ti,1}{i,iii}(ii,1)=mean(ask_vol_price_p{ti,1}{i,ii}(:,iii));
                    ask_vol_price_plus{ti,1}{i,iii}(ii,2)=size(ask_vol_price_p{ti,1}{i,ii}(:,iii),1);
                else
                    ask_vol_price_plus{ti,1}{i,iii}(ii,1)=0;
                    ask_vol_price_plus{ti,1}{i,iii}(ii,2)=0;
                end
                
                % bid_gap
                if isempty(bid_gap_m{ti,1}{i,ii})==0
                    bid_gap_minus{ti,1}{i,iii}(ii,1)=mean(bid_gap_m{ti,1}{i,ii}(:,iii));
                    bid_gap_minus{ti,1}{i,iii}(ii,2)=size(bid_gap_m{ti,1}{i,ii}(:,iii),1);
                else
                    bid_gap_minus{ti,1}{i,iii}(ii,1)=0;
                    bid_gap_minus{ti,1}{i,iii}(ii,2)=0;
                end
                if isempty(bid_gap_p{ti,1}{i,ii})==0
                    bid_gap_plus{ti,1}{i,iii}(ii,1)=mean(bid_gap_p{ti,1}{i,ii}(:,iii));
                    bid_gap_plus{ti,1}{i,iii}(ii,2)=size(bid_gap_p{ti,1}{i,ii}(:,iii),1);
                else
                    bid_gap_plus{ti,1}{i,iii}(ii,1)=0;
                    bid_gap_plus{ti,1}{i,iii}(ii,2)=0;
                end

                % ask_gap
                if isempty(ask_gap_m{ti,1}{i,ii})==0
                    ask_gap_minus{ti,1}{i,iii}(ii,1)=mean(ask_gap_m{ti,1}{i,ii}(:,iii));
                    ask_gap_minus{ti,1}{i,iii}(ii,2)=size(ask_gap_m{ti,1}{i,ii}(:,iii),1);
                else
                    ask_gap_minus{ti,1}{i,iii}(ii,1)=0;
                    ask_gap_minus{ti,1}{i,iii}(ii,2)=0;
                end
                if isempty(ask_gap_p{ti,1}{i,ii})==0
                    ask_gap_plus{ti,1}{i,iii}(ii,1)=mean(ask_gap_p{ti,1}{i,ii}(:,iii));
                    ask_gap_plus{ti,1}{i,iii}(ii,2)=size(ask_gap_p{ti,1}{i,ii}(:,iii),1);
                else
                    ask_gap_plus{ti,1}{i,iii}(ii,1)=0;
                    ask_gap_plus{ti,1}{i,iii}(ii,2)=0;
                end
            end
        end
    end
end


