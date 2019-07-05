function [] = lmc_flow_event_ver6_ABaover (firm_start,day_start,day_end,t)
% base code (2017.10.07)
% modified (2019.05.06)
% 9 + 1 + 2

addpath('/home/minyoung/data_minyoung/Research/Matlab_function/');
addpath('/home/minyoung/data_minyoung/Research/Price_Impact/');
date_seq=load('/home/minyoung/data_minyoung/Research/Price_Impact/date/date.txt');

% firm_start=1; %% firm~=10,14
firm_end=firm_start;

aggregate_size=1;
gamma_crit=0;
NI=9 + 1 + 2;

eps=10^(-4);

ind_day=0;

pvol_thre = logspace(4,7+log10(2),20);

base='/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/LSE_data2/';

for firm=firm_start:firm_end % remove [firm==73] case
    com_num=load('/home/minyoung/data_minyoung/Research/Price_Impact/top_73.csv');
    for day=day_start:day_end
            sprintf('%d, %d',firm,day)
%             a=load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/LSE_data2/%d/qth_corr_modi_%d_%d_%d_%d.csv',com_num(firm),com_num(firm),date_seq(day),aggregate_size,gamma_crit));
            a=load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/LSE_data2/%d/qth_corr_moc_%d_%d_%d_%d.csv',com_num(firm),com_num(firm),date_seq(day),aggregate_size,gamma_crit));
            
            name_trade=sprintf('%d/Trade_pre_%d_%d.csv',com_num(firm),com_num(firm),date_seq(day));
            name_trade_ask=sprintf('%d/Trade_pre_%d_ask_%d.csv',com_num(firm),com_num(firm),date_seq(day));
            name_trade_bid=sprintf('%d/Trade_pre_%d_bid_%d.csv',com_num(firm),com_num(firm),date_seq(day));
            name_trade_ask_vol=sprintf('%d/Trade_pre_%d_ask_vol_%d.csv',com_num(firm),com_num(firm),date_seq(day));
            name_trade_bid_vol=sprintf('%d/Trade_pre_%d_bid_vol_%d.csv',com_num(firm),com_num(firm),date_seq(day));

            file_trade=fullfile(base,name_trade);
            file_trade_ask=fullfile(base,name_trade_ask);
            file_trade_bid=fullfile(base,name_trade_bid);
            file_trade_ask_vol=fullfile(base,name_trade_ask_vol);
            file_trade_bid_vol=fullfile(base,name_trade_bid_vol);

            fid = fopen(file_trade, 'r');
            trade_temp = textscan(fid,'%f,%f,%f,%f,%f,%s');
            fclose(fid);

            trade_type = zeros(size(trade_temp{1,1},1),1); % 1905
            for i=1:size(trade_temp{1,6},1)
                type_temp = strsplit(char(trade_temp{1,6}(i,1)),',');
                if type_temp{1,2} == '0'
                    trade_type(i,1) = 0;    
                elseif type_temp{1,2} == 'AT'
                    trade_type(i,1) = 1;    
                elseif type_temp{1,2} == 'UT'
                    trade_type(i,1) = 2;    
                else
                    trade_type(i,1) = 3;    
                end
            end
            
            trade = zeros(size(trade_temp{1,1},1),4); % 1905
            for i=1:4
                trade(:,i)=trade_temp{1,i};
            end
            clear trade_temp
            
%             for i=1:size(trade,1)
%                 if trade(i,1)==trade(i+1)    
%             end
            

            fid = fopen(file_trade_ask, 'r');
            trade_ask = textscan(fid, '%s', 'Delimiter', '\n');
            fclose(fid);
            trade_ask = cellfun(@str2num, trade_ask{:}, 'UniformOutput', false);

            fid = fopen(file_trade_ask_vol, 'r');
            trade_ask_vol = textscan(fid, '%s', 'Delimiter', '\n');
            fclose(fid);
            trade_ask_vol = cellfun(@str2num, trade_ask_vol{:}, 'UniformOutput', false);

            fid = fopen(file_trade_bid, 'r');
            trade_bid = textscan(fid, '%s', 'Delimiter', '\n');
            fclose(fid);
            trade_bid = cellfun(@str2num, trade_bid{:}, 'UniformOutput', false);

            fid = fopen(file_trade_bid_vol, 'r');
            trade_bid_vol = textscan(fid, '%s', 'Delimiter', '\n');
            fclose(fid);
            trade_bid_vol = cellfun(@str2num, trade_bid_vol{:}, 'UniformOutput', false);

            clear file_trade file_trade_ask file_trade_ask_vol file_trade_bid file_trade_bid_vol
            clear name_trade name_trade_ask name_trade_ask_vol name_trade_bid name_trade_bid_vol
               
%             temp_ask = [];
%             for tt=1:length(trade_ask)
%                 for ttt=1:length(trade_ask{tt,1})
%                     temp_ask = [temp_ask; trade_ask{tt,1}(ttt,1)*ones(trade_ask_vol{tt,1}(ttt,1));];
%                 end
%             end
            
            trade_id=a(:,end);
            time_event=a(:,2);
            market=zeros(size(a,1),2);

            bid=a(:,end-2);
            ask=a(:,end-1);
            spread=ask-bid;
            price=(bid+ask)/2;

            market_temp=a(:,[3+2,4+2,9+2,10+2,15+2,16+2]);
            market(:,1)=sum(market_temp(:,[1,3,5]),2); % sell initiated
            market(:,2)=sum(market_temp(:,[2,4,6]),2); % buy initiated

            clear market_temp a

            price_return(1,1)=0;
            price_return(2:length(price),1)=price(2:end)-price(1:end-1);
            bid_return(1,1)=0;
            bid_return(2:length(price),1)=bid(2:end)-bid(1:end-1);
            ask_return(1,1)=0;
            ask_return(2:length(price),1)=ask(2:end)-ask(1:end-1);
            
            
            iter=1; % trading points except UT
            for i=1:length(trade(:,1))
                if trade(i,4)>0 && trade_type(i,1)~=2
                    trade_ind(iter,1)=i;
                    iter=iter+1;
                end
            end

            t_time=trade(trade_ind,1);
            t_time=unique(t_time);
%             t_time=t_time(2:end-1,:); % why 2:end-1 -> maybe removing UT trading

            % mark the trading event
            % overwrite trade_type at the same time
            for i=1:size(trade,1)
                for ii=1:size(t_time,1)
                    if trade(i,1)==t_time(ii,1)
                        trade(i,4)=1;
                        trade_type(i,1)=1;
                        break;
                    end
                end
            end
            trade(find(trade(:,4)>1+eps),4)=0; % delete the start & end point

            
            % delete the overlapping time (leave the last point -> leave the first point)
            clear trade_ind
            iter=1;
            for i=2:length(trade(:,1))
                if trade(i-1,1)~=trade(i,1)
                    trade_ind(iter,1)=i;
                    iter=iter+1;
                end
%                     if i==length(trade(:,1))
%                         trade_ind(iter,1)=i;
%                     end
            end


            trade=trade(trade_ind,:);
            trade_ask=trade_ask(trade_ind,1);
            trade_bid=trade_bid(trade_ind,1);
            trade_ask_vol=trade_ask_vol(trade_ind,1);
            trade_bid_vol=trade_bid_vol(trade_ind,1);
            trade_type=trade_type(trade_ind,1);
            clear trade_ind

            ind_day=ind_day+1;

            iter=1;
            for i=1:length(time_event)
                for ii=1:length(t_time)
                    if time_event(i,1)==t_time(ii)
                        time_ind(iter,1)=i;
                        iter=iter+1;
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            market = market(time_ind,:);
            price = price(time_ind,1);
            price_return = price_return(time_ind,1);
            bid_return = bid_return(time_ind,1);
            ask_return = ask_return(time_ind,1);
            trade_id = trade_id(time_ind,1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear time_ind

            
            bid_vol = zeros(size(trade,1),10);
            ask_vol = zeros(size(trade,1),10);
            bid_vol_price = zeros(size(trade,1),10);
            ask_vol_price = zeros(size(trade,1),10);
            bid_gap = zeros(size(trade,1),10);
            ask_gap = zeros(size(trade,1),10);
            bid_price = zeros(size(trade,1),1);
            ask_price = zeros(size(trade,1),1);
            % minimum order book size
            for i=1:size(trade,1)
                if size(trade_bid{i,1},2)>10 && size(trade_ask{i,1},2)>10
                    for ii=1:10
                        bid_vol(i,ii)=trade_bid_vol{i,1}(1,end-ii+1);
                        ask_vol(i,ii)=trade_ask_vol{i,1}(1,ii);
                        bid_vol_price(i,ii)=trade_bid{i,1}(1,end-ii+1);
                        ask_vol_price(i,ii)=trade_ask{i,1}(1,ii);
                        bid_gap(i,ii)=trade_bid{i,1}(1,end-ii+1)-trade_bid{i,1}(1,end-ii);
                        ask_gap(i,ii)=trade_ask{i,1}(1,ii+1)-trade_ask{i,1}(1,ii);
                    end
                else
                    for ii=1:10
                        bid_vol(i,ii)=0;
                        ask_vol(i,ii)=0;
                        bid_vol_price(i,ii)=0;
                        ask_vol_price(i,ii)=0;
                        bid_gap(i,ii)=0;
                        ask_gap(i,ii)=0;
                    end
                end
                bid_price(i,1)=trade_bid{i,1}(end);
                ask_price(i,1)=trade_ask{i,1}(1);
            end
            spread_a=ask_price-bid_price;
            clear trade_ask trade_bid trade_ask_vol trade_bid_vol

            % numbering the trading event
            iter=1;
            for i=1:size(trade,1)
                if trade(i,4)==1
                    trade(i,4)=iter;
                    iter=iter+1;
                end
            end       

            % [spread_b] is the spread at the event
            spread_b=zeros(size(market,1),1); % spread_b=zeros(size(limit,1),1); 1905
            for i=1:size(trade,1)
                if trade(i,4)>0
                    spread_b(trade(i,4),1)=trade(i,3)-trade(i,2);
                end
            end
            % [*_*_event] is *_* at the event
            ind=find(trade(:,4)>0);
            bid_vol_event=bid_vol(ind,:);
            ask_vol_event=ask_vol(ind,:);
            bid_vol_price_event=bid_vol_price(ind,:);
            ask_vol_price_event=ask_vol_price(ind,:);
            bid_gap_event=bid_gap(ind,:);
            ask_gap_event=ask_gap(ind,:);
            bid_price_event=bid_price(ind,1);
            ask_price_event=ask_price(ind,1);
            trade_type_event=trade_type(ind,1);
            
            
            % calculate the tick size
            bid_di=di(bid_price_event);
            ask_di=di(ask_price_event);
            tick_crit=[10^5,10000,5000,2000,1000,500,200,100,50,30,20,10,5,1,10^(-10)];
            for i=1:size(tick_crit,2)-1
                tick_bid_ind{i,1}=find(bid_price_event<tick_crit(i) & bid_price_event>=tick_crit(i+1));
                tick_ask_ind{i,1}=find(ask_price_event<tick_crit(i) & ask_price_event>=tick_crit(i+1));
                tick_di{i,1}=[abs(di(bid_price_event(find(bid_price_event<tick_crit(i) & bid_price_event>=tick_crit(i+1)))));abs(di(ask_price_event(find(ask_price_event<tick_crit(i) & ask_price_event>=tick_crit(i+1)))))];
            end
            for i=1:size(tick_crit,2)-1
                if isempty(min(tick_di{i,1}(find(tick_di{i,1}>0))))
                    tick_bin(i,1)=0;
                else
                    tick_bin(i,1)=min(tick_di{i,1}(find(tick_di{i,1}>0)));
                end
            end
            tick_bin_uni=unique(tick_bin);
            tick_bin_uni=tick_bin_uni(find(tick_bin_uni>0));
            tick_bin_uni_ind=[];
            for i=1:length(tick_bin_uni)
                tick_bin_uni_ind=[tick_bin_uni_ind;find(tick_bin==tick_bin_uni(i))];
            end
            tick_bin_uni_ind=sort(tick_bin_uni_ind);

            tick_bid_ind_valid=cell(length(tick_bin_uni),1);
            tick_ask_ind_valid=cell(length(tick_bin_uni),1);
            iter=0;
            for i=1:length(tick_bin)-1
                if tick_bin(i+1,1)>0 && tick_bin(i+1,1)~=tick_bin(i,1)
                    iter=iter+1;
                    tick_bin_valid(iter,1)=tick_bin(i+1,1);
                    tick_bid_ind_valid{iter,1}=[tick_bid_ind_valid{iter,1};tick_bid_ind{i+1,1}];
                    tick_ask_ind_valid{iter,1}=[tick_ask_ind_valid{iter,1};tick_ask_ind{i+1,1}];
                elseif tick_bin(i+1,1)>0 && tick_bin(i+1,1)==tick_bin(i,1)
                    tick_bin_valid(iter,1)=tick_bin(i+1,1);
                    tick_bid_ind_valid{iter,1}=[tick_bid_ind_valid{iter,1};tick_bid_ind{i+1,1}];
                    tick_ask_ind_valid{iter,1}=[tick_ask_ind_valid{iter,1};tick_ask_ind{i+1,1}];
                end
            end

            tick_bid=zeros(length(bid_price_event),1);
            tick_ask=zeros(length(ask_price_event),1);
            for i=1:length(tick_bin_valid)
                tick_bid(tick_bid_ind_valid{i,1},1)=tick_bin_valid(i,1);
                tick_ask(tick_ask_ind_valid{i,1},1)=tick_bin_valid(i,1);
            end

            % calculating the case of all spread
            clear ind_minus ind_plus
            ind_minus=cell(NI,1);
            ind_plus=cell(NI,1);
            for i=1:NI
                ind_minus{i,1}=cell(length(tick_bin_valid),1);
                ind_plus{i,1}=cell(length(tick_bin_valid),1);
            end
            
            market(find(t_time < 8.5*3600*1000),:) = 0;
            
            return_tick_ind = find(trade(:,4)>0);
            price_tick = (trade(:,2) + trade(:,3))/2;
            bid_tick = trade(:,2);
            ask_tick = trade(:,3);
            return_tick = price_tick(return_tick_ind+1) - price_tick(return_tick_ind);
            return_bid_tick = bid_tick(return_tick_ind+1) - bid_tick(return_tick_ind);
            return_ask_tick = ask_tick(return_tick_ind+1) - ask_tick(return_tick_ind);
            return_tick(find(t_time < 8.5*3600*1000),:) = 0;
            return_bid_tick(find(t_time < 8.5*3600*1000),:) = 0;
            return_ask_tick(find(t_time < 8.5*3600*1000),:) = 0;
%             return_tick_nonzero = return_tick(find(abs(return_tick)>0));
%             return_thre(1) = prctile(abs(return_tick_nonzero),50);
%             return_thre(2) = prctile(abs(return_tick_nonzero),70);
%             return_thre(3) = prctile(abs(return_tick_nonzero),90);
%             return_thre(4) = prctile(abs(return_tick_nonzero),99);

            
            % all spread
            % 1-9
            for i=1:length(tick_bin_valid)
                temp1=bid_vol_event(:,1);
                temp2=bid_vol_event(:,1);
                temp1=temp1;
                temp2=temp2+bid_vol_event(:,2);
                temp_onejump=intersect(find(temp1<=market(:,1) & temp2>market(:,1) & market(:,1)~=0),tick_bid_ind_valid{i,1});
                temp_twojump=intersect(find(temp2<=market(:,1) & market(:,1)~=0),tick_bid_ind_valid{i,1});
                                
                temp_onetick=intersect(find(return_bid_tick <= -1*tick_bin_valid(i) & return_bid_tick > -2*tick_bin_valid(i)),tick_bid_ind_valid{i,1});
                temp_twotick=intersect(find(return_bid_tick <= -2*tick_bin_valid(i) & return_bid_tick > -3*tick_bin_valid(i)),tick_bid_ind_valid{i,1});
                temp_threetick=intersect(find(return_bid_tick <= -3*tick_bin_valid(i) & return_bid_tick > -4*tick_bin_valid(i)),tick_bid_ind_valid{i,1});
                temp_fourtick=intersect(find(return_bid_tick <= -4*tick_bin_valid(i) & return_bid_tick > -5*tick_bin_valid(i)),tick_bid_ind_valid{i,1});
                temp_fivetick=intersect(find(return_bid_tick <= -5*tick_bin_valid(i) & return_bid_tick > -6*tick_bin_valid(i)),tick_bid_ind_valid{i,1});
                
                ind_minus{1,1}{i,1} = intersect(temp_onetick,temp_onejump);
                
                ind_minus{2,1}{i,1} = intersect(temp_twotick,temp_onejump);
                ind_minus{3,1}{i,1} = intersect(temp_twotick,temp_twojump);
                
                ind_minus{4,1}{i,1} = intersect(temp_threetick,temp_onejump);
                ind_minus{5,1}{i,1} = intersect(temp_threetick,temp_twojump);
                
                ind_minus{6,1}{i,1} = intersect(temp_fourtick,temp_onejump);
                ind_minus{7,1}{i,1} = intersect(temp_fourtick,temp_twojump);
                
                ind_minus{8,1}{i,1} = intersect(temp_fivetick,temp_onejump);
                ind_minus{9,1}{i,1} = intersect(temp_fivetick,temp_twojump);
             
                clear temp_onejump temp_twojump temp_threejump temp_fourjump temp_fivejump
                clear temp_onetick temp_twotick temp_threetick temp_fourtick temp_fivetick
                
                
                
                
                temp1=ask_vol_event(:,1);
                temp2=ask_vol_event(:,1);
                temp1=temp1;
                temp2=temp2+ask_vol_event(:,2);
                temp_onejump=intersect(find(temp1<=market(:,2) & temp2>market(:,2) & market(:,2)~=0),tick_bid_ind_valid{i,1});
                temp_twojump=intersect(find(temp2<=market(:,2) & market(:,2)~=0),tick_bid_ind_valid{i,1});
                  
                temp_onetick=intersect(find(return_ask_tick >= 1*tick_bin_valid(i) & return_ask_tick < 2*tick_bin_valid(i)),tick_bid_ind_valid{i,1});
                temp_twotick=intersect(find(return_ask_tick >= 2*tick_bin_valid(i) & return_ask_tick < 3*tick_bin_valid(i)),tick_bid_ind_valid{i,1});
                temp_threetick=intersect(find(return_ask_tick >= 3*tick_bin_valid(i) & return_ask_tick < 4*tick_bin_valid(i)),tick_bid_ind_valid{i,1});
                temp_fourtick=intersect(find(return_ask_tick >= 4*tick_bin_valid(i) & return_ask_tick < 5*tick_bin_valid(i)),tick_bid_ind_valid{i,1});
                temp_fivetick=intersect(find(return_ask_tick >= 5*tick_bin_valid(i) & return_ask_tick < 6*tick_bin_valid(i)),tick_bid_ind_valid{i,1});
                
                ind_plus{1,1}{i,1} = intersect(temp_onetick,temp_onejump);
                
                ind_plus{2,1}{i,1} = intersect(temp_twotick,temp_onejump);
                ind_plus{3,1}{i,1} = intersect(temp_twotick,temp_twojump);
                
                ind_plus{4,1}{i,1} = intersect(temp_threetick,temp_onejump);
                ind_plus{5,1}{i,1} = intersect(temp_threetick,temp_twojump);
                
                ind_plus{6,1}{i,1} = intersect(temp_fourtick,temp_onejump);
                ind_plus{7,1}{i,1} = intersect(temp_fourtick,temp_twojump);
                
                ind_plus{8,1}{i,1} = intersect(temp_fivetick,temp_onejump);
                ind_plus{9,1}{i,1} = intersect(temp_fivetick,temp_twojump);
             
                clear temp_onejump temp_twojump temp_threejump temp_fourjump temp_fivejump
                clear temp_onetick temp_twotick temp_threetick temp_fourtick temp_fivetick

            end     
            % 10
            for i=1:length(tick_bin_valid)
                temp1=bid_vol_event(:,1);
                temp2=bid_vol_event(:,1);
                temp1=temp1;
                temp2=temp2+bid_vol_event(:,2);
                temp_nojump=intersect(find(temp1>market(:,1) & market(:,1)~=0),tick_bid_ind_valid{i,1});            
                ind_minus{10,1}{i,1} = temp_nojump;
             
                clear temp_onejump temp_twojump temp_threejump temp_fourjump temp_fivejump
                clear temp_onetick temp_twotick temp_threetick temp_fourtick temp_fivetick
                clear temp_nojump            
                
                
                temp1=ask_vol_event(:,1);
                temp2=ask_vol_event(:,1);
                temp1=temp1;
                temp2=temp2+ask_vol_event(:,2);
                temp_nojump=intersect(find(temp1>market(:,2) & market(:,2)~=0),tick_bid_ind_valid{i,1});
                ind_plus{10,1}{i,1} = temp_nojump;
                clear temp_onejump temp_twojump temp_threejump temp_fourjump temp_fivejump
                clear temp_onetick temp_twotick temp_threetick temp_fourtick temp_fivetick
                clear temp_nojump
            end     
            
            % 11-12
            for i=1:length(tick_bin_valid)
                temp1=bid_vol_event(:,1);
                temp2=bid_vol_event(:,1);
                temp1=temp1;
                temp2=temp2+bid_vol_event(:,2);
                temp_onejump=intersect(find(temp1<=market(:,1) & temp2>market(:,1) & market(:,1)~=0),tick_bid_ind_valid{i,1});
                temp_twojump=intersect(find(temp2<=market(:,1) & market(:,1)~=0),tick_bid_ind_valid{i,1});
                                
                temp_onetick=intersect(find(return_bid_tick <= -1*tick_bin_valid(i) & return_bid_tick > -2*tick_bin_valid(i)),tick_bid_ind_valid{i,1});
                temp_twotick=intersect(find(return_bid_tick <= -2*tick_bin_valid(i) & return_bid_tick > -3*tick_bin_valid(i)),tick_bid_ind_valid{i,1});
                temp_threetick=intersect(find(return_bid_tick <= -3*tick_bin_valid(i) & return_bid_tick > -4*tick_bin_valid(i)),tick_bid_ind_valid{i,1});
                temp_fourtick=intersect(find(return_bid_tick <= -4*tick_bin_valid(i) & return_bid_tick > -5*tick_bin_valid(i)),tick_bid_ind_valid{i,1});
                temp_fivetick=intersect(find(return_bid_tick <= -5*tick_bin_valid(i) & return_bid_tick > -6*tick_bin_valid(i)),tick_bid_ind_valid{i,1});
                temp_overtick=intersect(find(return_bid_tick <= -6*tick_bin_valid(i)),tick_bid_ind_valid{i,1});

                
                ind_minus{11,1}{i,1} = intersect(temp_overtick,temp_onejump);
                ind_minus{12,1}{i,1} = intersect(temp_overtick,temp_twojump);
                clear temp_onejump temp_twojump temp_threejump temp_fourjump temp_fivejump
                clear temp_onetick temp_twotick temp_threetick temp_fourtick temp_fivetick
                clear temp_overtick
                
                
                
                
                temp1=ask_vol_event(:,1);
                temp2=ask_vol_event(:,1);
                temp1=temp1;
                temp2=temp2+ask_vol_event(:,2);
                temp_onejump=intersect(find(temp1<=market(:,2) & temp2>market(:,2) & market(:,2)~=0),tick_bid_ind_valid{i,1});
                temp_twojump=intersect(find(temp2<=market(:,2) & market(:,2)~=0),tick_bid_ind_valid{i,1});
                  
                temp_onetick=intersect(find(return_ask_tick >= 1*tick_bin_valid(i) & return_ask_tick < 2*tick_bin_valid(i)),tick_bid_ind_valid{i,1});
                temp_twotick=intersect(find(return_ask_tick >= 2*tick_bin_valid(i) & return_ask_tick < 3*tick_bin_valid(i)),tick_bid_ind_valid{i,1});
                temp_threetick=intersect(find(return_ask_tick >= 3*tick_bin_valid(i) & return_ask_tick < 4*tick_bin_valid(i)),tick_bid_ind_valid{i,1});
                temp_fourtick=intersect(find(return_ask_tick >= 4*tick_bin_valid(i) & return_ask_tick < 5*tick_bin_valid(i)),tick_bid_ind_valid{i,1});
                temp_fivetick=intersect(find(return_ask_tick >= 5*tick_bin_valid(i) & return_ask_tick < 6*tick_bin_valid(i)),tick_bid_ind_valid{i,1});
                temp_overtick=intersect(find(return_ask_tick >= 6*tick_bin_valid(i)),tick_bid_ind_valid{i,1});

                
                ind_plus{11,1}{i,1} = intersect(temp_overtick,temp_onejump);
                ind_plus{12,1}{i,1} = intersect(temp_overtick,temp_twojump);
                clear temp_onejump temp_twojump temp_threejump temp_fourjump temp_fivejump
                clear temp_onetick temp_twotick temp_threetick temp_fourtick temp_fivetick
                clear temp_overtick
            end     
            
            ind_minus_temp = ind_minus;
            ind_plus_temp = ind_plus;
            clear ind_minus ind_plus
            
            pv(:,1) = bid_price_event .* market(:,1); % sell initiated
            pv(:,2) = ask_price_event .* market(:,2); % buy initiated
            
%             clear pv pv_thre_1 pv_thre_2
            
%             for large_ind = 0:1
            for pv_ind = 1:length(pv_thre)

                pv_thre_1 = find(pv(:,1) < pvol_thre(pv_ind) & pv(:,1) > 0);
                pv_thre_2 = find(pv(:,2) < pvol_thre(pv_ind) & pv(:,2) > 0);

                for i = 1:NI
                    for ii = 1:length(ind_minus_temp{i,1})
                        ind_minus{i,1}{ii,1} = intersect(ind_minus_temp{i,1}{ii,1}, pv_thre_1); % sell initiated
                        ind_plus{i,1}{ii,1} = intersect(ind_plus_temp{i,1}{ii,1}, pv_thre_2); % buy initiated
                    end
                end
                clear pv_thre_1 pv_thre_2


                x_time_bin = linspace(8*3600*1000,16.5*3600*1000,510);
                for i=1:NI
                    for ii=1:length(ind_minus{1,1})
                        [k1{i,ii} emp] = hist(t_time(eval(sprintf('ind_minus{%d,1}{%d,1}',i,ii))),x_time_bin);
                    end
                end
                clear emp
                for i=1:NI
                    for ii=1:length(ind_minus{1,1})
                        [k2{i,ii} emp] = hist(t_time(eval(sprintf('ind_plus{%d,1}{%d,1}',i,ii))),x_time_bin);
                    end
                end
                clear emp

                [spread_plus,spread_minus,bid_vol_minus,bid_vol_plus,ask_vol_minus,ask_vol_plus,bid_vol_price_minus,bid_vol_price_plus,ask_vol_price_minus,ask_vol_price_plus,bid_gap_minus,bid_gap_plus,ask_gap_minus,ask_gap_plus,bid_price_minus,bid_price_plus,ask_price_minus,ask_price_plus]=calculating_event_ver6_1905(t,trade,spread_a,bid_vol,ask_vol,bid_vol_price,ask_vol_price,bid_gap,ask_gap,bid_price,ask_price,tick_bin_valid,ind_minus,ind_plus,NI);
                if ~exist(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result_AB_1905/%d',com_num(firm)),'dir')
                    mkdir(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result_AB_1905/%d',com_num(firm)))
                end

                save(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result_AB_1905/%d/lmc_flow_event_ver6_AB_k_1905_%d_%d_%d_%d_%d.mat',com_num(firm),com_num(firm),day,1,t,pv_ind),'k1','k2')
                save(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result_AB_1905/%d/lmc_flow_event_ver6_AB_sign_1905_%d_%d_%d_%d_%d.mat',com_num(firm),com_num(firm),day,1,t,pv_ind),'spread_plus','spread_minus','bid_vol_minus','bid_vol_plus','ask_vol_minus','ask_vol_plus','bid_vol_price_minus','bid_vol_price_plus','ask_vol_price_minus','ask_vol_price_plus','bid_gap_minus','bid_gap_plus','ask_gap_minus','ask_gap_plus','bid_price_minus','bid_price_plus','ask_price_minus','ask_price_plus','tick_bin_valid')
                save(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result_AB_1905/%d/lmc_flow_event_ver6_AB_sign_ind_1905_%d_%d_%d_%d_%d.mat',com_num(firm),com_num(firm),day,1,t,pv_ind),'ind_minus','ind_plus')
                save(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result_AB_1905/%d/lmc_flow_event_ver6_AB_sign_var_1905_%d_%d_%d_%d_%d.mat',com_num(firm),com_num(firm),day,1,t,pv_ind),'ask_gap_event','ask_price_event','ask_return','ask_vol_event','bid_gap_event','bid_price_event','bid_return','bid_vol_event','market','price','price_return','t_time','trade_id')

                clear k1 k2
                clear ind_minus ind_plus
                clear spread_plus spread_minus bid_vol_minus bid_vol_plus ask_vol_minus 
                clear ask_vol_plus bid_vol_price_minus bid_vol_price_plus ask_vol_price_minus ask_vol_price_plus
                clear bid_gap_minus bid_gap_plus ask_gap_minus ask_gap_plus bid_price_minus bid_price_plus 
                clear ask_price_minus ask_price_plus

            end
%             end
            clear pv
            clear ind_minus_temp ind_plus_temp

            
            clear bid_di ask_di tick_di
            clear ind_bid ind_ask
            clear temp1 temp2
            clear bid_vol_price_minus bid_vol_price_plus ask_vol_price_minus ask_vol_price_plus
            clear spread_plus spread_minus bid_vol_minus bid_vol_plus ask_vol_minus ask_vol_plus bid_gap_minus bid_gap_plus ask_gap_minus
            clear ask_gap_plus bid_price_minus bid_price_plus ask_price_minus ask_price_plus flow_minus flow_plus update_time_gap_minus update_time_gap_plus
            clear bid_vol_price_m_avg bid_vol_price_p_avg ask_vol_price_m_avg ask_vol_price_p_avg
            clear bid_di ask_di tick_di
            clear ind_bid ind_ask
            clear temp1 temp2
            clear bid_vol_price_minus bid_vol_price_plus ask_vol_price_minus ask_vol_price_plus
            clear spread_plus spread_minus bid_vol_minus bid_vol_plus ask_vol_minus ask_vol_plus bid_gap_minus bid_gap_plus ask_gap_minus
            clear ask_gap_plus bid_price_minus bid_price_plus ask_price_minus ask_price_plus flow_minus flow_plus update_time_gap_minus update_time_gap_plus
            clear bid_vol_price_m_avg bid_vol_price_p_avg ask_vol_price_m_avg ask_vol_price_p_avg
            clear flow
            clear limit market cancel lc_plus orderflow spread time bid_vol ask_vol bid_gap ask_gap
            clear bid_vol_price ask_vol_price bid_vol_price_event ask_vol_price_event
            clear ask bid cancel limit price_return market price spread
            clear limit market cancel lc_plus orderflow spread time bid_vol ask_vol bid_gap ask_gap
            clear ask bid cancel limit price_return market price spread x y
            clear trade trade_ask trade_bid trade_ask_vol trade_bid_vol bid_return ask_return
            clear ask_gap_event ask_vol_event ask_price bid_gap_event bid_vol_event bid_price ind spread_a spread_b t_time
            clear ask_price_event bid_price_event spread_crit_ask spread_crit_bid tick_bid tick_ask tick_bin_valid tick_ask_ind_valid tick_bid_ind_valid
            clear tick_ask_ind tick_bid_ind tick_bin
            clear trade_id lc_minus time_event time_gap trade_type trade_type_event type_temp
            clear update_time_gap
            clear price_tick return_tick return_thre
            clear ask_tick bid_tick return_ask_tick return_bid_tick return_tick_ind
    end
end





