function [] = burst_individual_save_1906(firm_start,firm_end,win_size,mv_size,burst_ind)
% tracking burstiness, memory coefficient

addpath('/home/minyoung/data_minyoung/Research/Matlab_function/');
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/randraw/');
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/arrow/');
addpath('/home/minyoung/data_minyoung/Research/Price_Impact/');
date_seq = load('/home/minyoung/data_minyoung/Research/Price_Impact/date/date.txt');
com_num = load('/home/minyoung/data_minyoung/Research/Price_Impact/top_69.csv');
%com_num = load('/home/minyoung/data_minyoung/Research/Price_Impact/top_73.csv');

% win_size = 10;
% mv_size = 1;
id = 1;
t = 1000;
NI = 4;
% firm_end = firm_start;
len_firm = firm_end - firm_start + 1;

mem_coef = cell(ceil((169-win_size+1)/mv_size),2);
burst = cell(ceil((169-win_size+1)/mv_size),2);
for i=1:length(mem_coef)
    mem_coef{i,1} = zeros(len_firm,NI);
    burst{i,1} = zeros(len_firm,NI);
    mem_coef{i,2} = zeros(len_firm,NI);
    burst{i,2} = zeros(len_firm,NI);
end

volatility = zeros(ceil((169-win_size+1)/mv_size),len_firm);
trading_volume = zeros(ceil((169-win_size+1)/mv_size),len_firm);
price_df = zeros(ceil((169-win_size+1)/mv_size),len_firm); 

for day_loop=1:ceil((169-win_size+1)/mv_size)
    day_loop

    day_start = 1 + (day_loop-1)*mv_size;
    day_end = day_start + win_size -1; 
    
    if win_size == 30 && mv_size == 30 && day_loop == ceil((169-win_size+1)/mv_size)
        day_start = 140;
        day_end = 169;
    end
    
    if day_start<=26 && day_end>=26 % 26th(20080908),new trading flatform with MS
        len_day = day_end - day_start + 1 -1; % remove 26th
    else
        len_day = day_end - day_start + 1;
    end
    if len_day ==0
        len_day = 1;    
    end
    
    interevent_time_minus = cell(len_day,len_firm);
    interevent_time_plus = cell(len_day,len_firm);
    firm_ind = 1;
    for firm=firm_start:firm_end
        day_ind = 1;
        for day=day_start:day_end
            if day~=26
                interevent_time_minus{day_ind,firm_ind} = cell(NI,1);
                interevent_time_plus{day_ind,firm_ind} = cell(NI,1);
                day_ind = day_ind + 1;
            end
        end
        firm_ind = firm_ind + 1;
    end
    time_minus = cell(len_day,len_firm);
    time_plus = cell(len_day,len_firm);
    firm_ind = 1;
    for firm=firm_start:firm_end
        day_ind = 1;
        for day=day_start:day_end
            if day~=26
                time_minus{day_ind,firm_ind} = cell(NI,1);
                time_plus{day_ind,firm_ind} = cell(NI,1);
                day_ind = day_ind + 1;
            end
        end
        firm_ind = firm_ind + 1;
    end

    % need to aggregate length(ind_minus{1,1})
    firm_ind = 1;
    for firm=firm_start:firm_end
%         firm
        volume = [];
        abs_return = [];
        price = [];

        day_ind = 1;
        for day=day_start:day_end
            if day~=26
                load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result_15/lmc_flow_event_ver6_repair_g_time_%d_%d_%d_%d.mat',com_num(firm_ind),day,1,t))
                load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_AB_sign_ind_%d_%d_%d_%d.mat',com_num(firm_ind),day,id,t))
%                 load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_AB_sign_%d_%d_%d_%d.mat',com_num(firm_ind),day,id,t))
                AB_var = load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_AB_sign_var_%d_%d_%d_%d.mat',com_num(firm_ind),day,id,t),'price','market');
                ABa = load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_ABa_sign_ind_%d_%d_%d_%d.mat',com_num(firm_ind),day,id,t));
                AB_over = load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result_over/lmc_flow_event_ver6_AB_over_sign_ind_%d_%d_%d_%d.mat',com_num(firm_ind),day,id,t));
                price_return = log(AB_var.price(2:end,1)) - log(AB_var.price(1:end-1,1));
                price = [price; AB_var.price(end,1)];
                volume = [volume; AB_var.market];
                abs_return = [abs_return; abs(price_return)];
                if firm == 60 && day == 103 % correct error
                    abs_return = abs_return(1:end-1);    
                end
                clear AB_var price_return
            
                for ii=1:length(ind_minus{1,1})
                    ind_minus{10,1}{ii,1} = ABa.ind_minus{1,1}{ii,1};
                    ind_plus{10,1}{ii,1} = ABa.ind_plus{1,1}{ii,1};   
                    
                    ind_minus{11,1}{ii,1} = AB_over.ind_minus{1,1}{ii,1};
                    ind_plus{11,1}{ii,1} = AB_over.ind_plus{1,1}{ii,1};   
                    ind_minus{12,1}{ii,1} = AB_over.ind_minus{2,1}{ii,1};
                    ind_plus{12,1}{ii,1} = AB_over.ind_plus{2,1}{ii,1};   
                end

                for ii=1:length(ind_minus{1,1})    
% %                     ind_minus_temp{1,1}{ii,1} = sort(ind_minus{10,1}{ii,1});
                    ind_minus_temp{1,1}{ii,1} = sort([ind_minus{10,1}{ii,1};ind_plus{10,1}{ii,1}]);
                    %ind_minus_temp{2,1}{ii,1} = sort([ind_minus{1,1}{ii,1};ind_minus{2,1}{ii,1};ind_minus{4,1}{ii,1};ind_minus{6,1}{ii,1};ind_minus{8,1}{ii,1};ind_minus{11,1}{ii,1}]);
% %                     ind_minus_temp{2,1}{ii,1} = sort([ind_minus{2,1}{ii,1};ind_minus{4,1}{ii,1};ind_minus{6,1}{ii,1};ind_minus{8,1}{ii,1};ind_minus{11,1}{ii,1}]);
% %                     ind_minus_temp{3,1}{ii,1} = sort([ind_minus{3,1}{ii,1};ind_minus{5,1}{ii,1};ind_minus{7,1}{ii,1};ind_minus{9,1}{ii,1};ind_minus{12,1}{ii,1}]);
                    ind_minus_temp{2,1}{ii,1} = sort([ind_minus{2,1}{ii,1};ind_minus{4,1}{ii,1};ind_minus{6,1}{ii,1};ind_minus{8,1}{ii,1};ind_minus{11,1}{ii,1};ind_plus{2,1}{ii,1};ind_plus{4,1}{ii,1};ind_plus{6,1}{ii,1};ind_plus{8,1}{ii,1};ind_plus{11,1}{ii,1}]);
                    ind_minus_temp{3,1}{ii,1} = sort([ind_minus{3,1}{ii,1};ind_minus{5,1}{ii,1};ind_minus{7,1}{ii,1};ind_minus{9,1}{ii,1};ind_minus{12,1}{ii,1};ind_plus{3,1}{ii,1};ind_plus{5,1}{ii,1};ind_plus{7,1}{ii,1};ind_plus{9,1}{ii,1};ind_plus{12,1}{ii,1}]);
                    ind_minus_temp{4,1}{ii,1} = sort([ind_minus{1,1}{ii,1};ind_plus{1,1}{ii,1}]);
%                     ind_minus_temp{4,1}{ii,1} = sort([ind_minus_temp{1,1}{ii,1};ind_minus_temp{2,1}{ii,1};ind_minus_temp{3,1}{ii,1}]); %%%%%%%%%%%%%%%%%

                    ind_plus_temp{1,1}{ii,1} = sort(ind_plus{10,1}{ii,1});
                    %ind_plus_temp{2,1}{ii,1} = sort([ind_plus{1,1}{ii,1};ind_plus{2,1}{ii,1};ind_plus{4,1}{ii,1};ind_plus{6,1}{ii,1};ind_plus{8,1}{ii,1};ind_plus{11,1}{ii,1}]);
                    ind_plus_temp{2,1}{ii,1} = sort([ind_plus{2,1}{ii,1};ind_plus{4,1}{ii,1};ind_plus{6,1}{ii,1};ind_plus{8,1}{ii,1};ind_plus{11,1}{ii,1}]);
                    ind_plus_temp{3,1}{ii,1} = sort([ind_plus{3,1}{ii,1};ind_plus{5,1}{ii,1};ind_plus{7,1}{ii,1};ind_plus{9,1}{ii,1};ind_plus{12,1}{ii,1}]);
                    ind_plus_temp{4,1}{ii,1} = sort(ind_plus{1,1}{ii,1});
%                     ind_plus_temp{4,1}{ii,1} = sort([ind_plus_temp{1,1}{ii,1};ind_plus_temp{2,1}{ii,1};ind_plus_temp{3,1}{ii,1}]);
                    
%                     ind_minus_temp{5,1}{ii,1} = sort(unique([ind_minus_temp{4,1}{ii,1};ind_plus_temp{4,1}{ii,1}]));
%                     ind_plus_temp{5,1}{ii,1} = sort(unique([ind_minus_temp{4,1}{ii,1};ind_plus_temp{4,1}{ii,1}]));                    
                end
                clear ind_minus ind_plus
                ind_minus = ind_minus_temp;
                ind_plus = ind_plus_temp;
                clear ind_minus_temp ind_plus_temp
            

                cut_t = 1;
                start_time = (8+cut_t)*3600*1000;
                end_time = (16.5-cut_t)*3600*1000;
                for i=1:NI
                    for ii=1:length(ind_minus{1,1})
                        ind_series1 = find(time_event>start_time & time_event<end_time);
                        ind_series2 = find(ind_minus{i,1}{ii,1}>min(ind_series1) & ind_minus{i,1}{ii,1}<max(ind_series1));
                        ind_minus{i,1}{ii,1} = ind_minus{i,1}{ii,1}(ind_series2);
                        clear ind_series1 ind_series2
                        ind_series1 = find(time_event>start_time & time_event<end_time);
                        ind_series2 = find(ind_plus{i,1}{ii,1}>min(ind_series1) & ind_plus{i,1}{ii,1}<max(ind_series1));
                        ind_plus{i,1}{ii,1} = ind_plus{i,1}{ii,1}(ind_series2);
                        clear ind_series1 ind_series2
                    end
                end


                for i=1:NI
                    for ii=1:length(ind_minus{1,1})
                        if ii==1
                            interevent_time_minus{day_ind,firm_ind}{i,1} = cell(length(ind_minus{1,1}),1);
                            interevent_time_plus{day_ind,firm_ind}{i,1} = cell(length(ind_plus{1,1}),1);
                            time_minus{day_ind,firm_ind}{i,1} = cell(length(ind_minus{1,1}),1);
                            time_plus{day_ind,firm_ind}{i,1} = cell(length(ind_plus{1,1}),1);
                        end
                        if length(ind_minus{i,1}{ii,1})>1 && length(ind_plus{i,1}{ii,1})>1
                            interevent_time_minus{day_ind,firm_ind}{i,1}{ii,1} = time_event(ind_minus{i,1}{ii,1}(2:end)) - time_event(ind_minus{i,1}{ii,1}(1:end-1));
                            interevent_time_plus{day_ind,firm_ind}{i,1}{ii,1} = time_event(ind_plus{i,1}{ii,1}(2:end)) - time_event(ind_plus{i,1}{ii,1}(1:end-1));
                        else
                            interevent_time_minus{day_ind,firm_ind}{i,1}{ii,1} = [];
                            interevent_time_plus{day_ind,firm_ind}{i,1}{ii,1} = [];
                        end
                        time_minus{day_ind,firm_ind}{i,1}{ii,1} = time_event(ind_minus{i,1}{ii,1});
                        time_plus{day_ind,firm_ind}{i,1}{ii,1} = time_event(ind_plus{i,1}{ii,1});
                    end
                end
                clear ind_minus ind_plus time_event
                day_ind = day_ind + 1;
            end
        end
        
        if isempty(abs_return) % if len=1, day_loop=26
            abs_return = 0;
        end
        if isempty(volume) % if len=1, day_loop=26
            volume = 0;
        end
        if isempty(price) % if len=1, day_loop=26
            price = 0;
        end
        volatility(day_loop,firm_ind) = std(abs_return);
        trading_volume(day_loop,firm_ind) = sum(sum(volume));
        price_df(day_loop,firm_ind) = mean(price);
        firm_ind = firm_ind + 1;
    end

    firm_ind = 1;
    for firm=firm_start:firm_end
        interevent_time_minus_ag = cell(NI,firm_ind);
        interevent_time_plus_ag = cell(NI,firm_ind);
        firm_ind = firm_ind + 1;
    end

    firm_ind = 1;
    for firm=firm_start:firm_end
        day_ind = 1;
        for day=day_start:day_end
            if day~=26
                for i=1:NI
                    for ii=1:length(interevent_time_minus{day_ind,ii}{1,1})
                        interevent_time_minus_ag{i,firm_ind} = [interevent_time_minus_ag{i,firm_ind}; interevent_time_minus{day_ind,firm_ind}{i,1}{ii,1}];
                        interevent_time_plus_ag{i,firm_ind} = [interevent_time_plus_ag{i,firm_ind}; interevent_time_plus{day_ind,firm_ind}{i,1}{ii,1}];
                    end
                end
                day_ind = day_ind + 1;
            end
        end
        firm_ind = firm_ind + 1;
    end
    clear interevent_time_minus interevent_time_plus

    % % normalization
    % for firm=firm_start:firm_end
    %     for i=1:NI
    %         interevent_time_minus_ag{i,firm} = (interevent_time_minus_ag{i,firm})/std(interevent_time_minus_ag{i,firm});
    %         interevent_time_plus_ag{i,firm} = (interevent_time_plus_ag{i,firm})/std(interevent_time_plus_ag{i,firm});
    %     end
    % end




    f_ind = 1;
    
    for f=1:len_firm
        for i=1:NI
            if length(interevent_time_minus_ag{i,f_ind})>5 && length(interevent_time_plus_ag{i,f_ind})>5
                
                % minus
%                 mean_iet = mean(interevent_time_minus_ag{i,f_ind});
%                 std_iet = std(interevent_time_minus_ag{i,f_ind});

                y = interevent_time_minus_ag{i,f_ind};
%                 y = interevent_time_minus_ag{i,f_ind}/mean_iet;
                mem_coef{day_loop,1}(f_ind,i) = memory_coef(y);
                burst{day_loop,1}(f_ind,i) = burstiness(y,burst_ind);
                
                % plus
%                 mean_iet = mean(interevent_time_plus_ag{i,f_ind});
%                 std_iet = std(interevent_time_plus_ag{i,f_ind});

                y = interevent_time_plus_ag{i,f_ind};
%                 y = interevent_time_minus_ag{i,f_ind}/mean_iet;
                mem_coef{day_loop,2}(f_ind,i) = memory_coef(y);
                burst{day_loop,2}(f_ind,i) = burstiness(y,burst_ind);
            end
        end
        f_ind = f_ind + 1;
    end
    clear interevent_time_minus_ag interevent_time_plus_ag
    clear time_minus time_plus
end


for i=1:size(burst,1)
    if win_size == 1
        if i==26
            for ii=1:NI
                mean_burst_minus(i,ii) = 0;
                mean_mem_minus(i,ii) = 0;

                mean_burst_plus(i,ii) = 0;
                mean_mem_plus(i,ii) = 0;

                std_burst_minus(i,ii) = 0;
                std_mem_minus(i,ii) = 0;

                std_burst_plus(i,ii) = 0;
                std_mem_plus(i,ii) = 0;
            end    
        else
            for ii=1:NI
                mean_burst_minus(i,ii) = mean(burst{i,1}(:,ii));
                mean_mem_minus(i,ii) = mean(mem_coef{i,1}(:,ii));

                mean_burst_plus(i,ii) = mean(burst{i,2}(:,ii));
                mean_mem_plus(i,ii) = mean(mem_coef{i,2}(:,ii));

                std_burst_minus(i,ii) = std(burst{i,1}(:,ii));
                std_mem_minus(i,ii) = std(mem_coef{i,1}(:,ii));

                std_burst_plus(i,ii) = std(burst{i,2}(:,ii));
                std_mem_plus(i,ii) = std(mem_coef{i,2}(:,ii));
            end
        end
    elseif win_size ~= 1
        for ii=1:NI
            mean_burst_minus(i,ii) = mean(burst{i,1}(:,ii));
            mean_mem_minus(i,ii) = mean(mem_coef{i,1}(:,ii));

            mean_burst_plus(i,ii) = mean(burst{i,2}(:,ii));
            mean_mem_plus(i,ii) = mean(mem_coef{i,2}(:,ii));

            std_burst_minus(i,ii) = std(burst{i,1}(:,ii));
            std_mem_minus(i,ii) = std(mem_coef{i,1}(:,ii));

            std_burst_plus(i,ii) = std(burst{i,2}(:,ii));
            std_mem_plus(i,ii) = std(mem_coef{i,2}(:,ii));
        end
    end
end


ks_alpha = .01;
for i=1:size(burst,1)
    ks_burst_minus(i,1) = kstest2(burst{i,1}(:,1),burst{i,1}(:,2),'Alpha',ks_alpha);
    ks_burst_minus(i,2) = kstest2(burst{i,1}(:,1),burst{i,1}(:,3),'Alpha',ks_alpha);
    ks_burst_minus(i,3) = kstest2(burst{i,1}(:,2),burst{i,1}(:,3),'Alpha',ks_alpha);

    ks_burst_plus(i,1) = kstest2(burst{i,2}(:,1),burst{i,2}(:,2),'Alpha',ks_alpha);
    ks_burst_plus(i,2) = kstest2(burst{i,2}(:,1),burst{i,2}(:,3),'Alpha',ks_alpha);
    ks_burst_plus(i,3) = kstest2(burst{i,2}(:,2),burst{i,2}(:,3),'Alpha',ks_alpha);
end
for i=1:size(mem_coef,1)
    ks_mem_minus(i,1) = kstest2(mem_coef{i,1}(:,1),mem_coef{i,1}(:,2),'Alpha',ks_alpha);
    ks_mem_minus(i,2) = kstest2(mem_coef{i,1}(:,1),mem_coef{i,1}(:,3),'Alpha',ks_alpha);
    ks_mem_minus(i,3) = kstest2(mem_coef{i,1}(:,2),mem_coef{i,1}(:,3),'Alpha',ks_alpha);

    ks_mem_plus(i,1) = kstest2(mem_coef{i,2}(:,1),mem_coef{i,2}(:,2),'Alpha',ks_alpha);
    ks_mem_plus(i,2) = kstest2(mem_coef{i,2}(:,1),mem_coef{i,2}(:,3),'Alpha',ks_alpha);
    ks_mem_plus(i,3) = kstest2(mem_coef{i,2}(:,2),mem_coef{i,2}(:,3),'Alpha',ks_alpha);
end

mean_volatility = mean(volatility,2);
mean_trading_volume = mean(trading_volume,2);
mean_price = mean(price_df,2);


save(sprintf('./burst_cal/burst_individual_1906_%d_%d_%d.mat',win_size,mv_size,burst_ind),'burst','ks_burst_minus','ks_burst_plus','ks_mem_minus','ks_mem_plus','mean_burst_minus','mean_burst_plus','mean_mem_minus','mean_mem_plus','mean_trading_volume','mean_volatility','mean_price','mem_coef','std_burst_minus','std_burst_plus','std_mem_minus','std_mem_plus','trading_volume','volatility')

%%
% 
% for i = 1:160
%     for ii = 1:69
%         for iii = 1:3
%             for iv = 1:2
%                 b{ii,iv}(i,iii) = burst{i,iv}(ii,iii);
%                 m{ii,iv}(i,iii) = mem_coef{i,iv}(ii,iii);
%             end
%         end
%     end
% end
% 
% %%
% for li = 1:4
%     ind = 1;
%     figure
%     for i = (li-1)*24 + 1: li*24
%         subplot(6,4,ind)
%         plot(b{i,1}(:,3))
%         ind = ind + 1;
%     end
% end
% 
% %%
% burst_timeseries = zeros(160,3);
% mem_coef_timeseries = zeros(160,3);
% for i = 1:160
%     for ii = 1:3
%         burst_timeseries(i,ii) = mean(burst{i,1}(:,ii));
%         mem_coef_timeseries(i,ii) = mean(mem_coef{i,1}(:,ii));
%         
%         std_burst_timeseries(i,ii) = std(burst{i,1}(:,ii));
%         std_mem_coef_timeseries(i,ii) = std(mem_coef{i,1}(:,ii));
%         
% %         burst_timeseries(i,ii) = mean(delniz(burst{i,1}(:,ii),3,1),1);
% %         mem_coef_timeseries(i,ii) = mean(delniz(mem_coef{i,1}(:,ii),3,1),1);
% %         
% %         std_burst_timeseries(i,ii) = std(delniz(burst{i,1}(:,ii),3,1),1);
% %         std_mem_coef_timeseries(i,ii) = std(delniz(mem_coef{i,1}(:,ii),3,1),1);
%     end
% end
% 
% figure
% subplot(1,2,1)
% plot(burst_timeseries)
% subplot(1,2,2)
% plot(mem_coef_timeseries)
% 
% 
% cif = tinv(0.975, 69-1) / sqrt(69);
% mean_burst_minus = burst_timeseries;
% mean_mem_minus = mem_coef_timeseries;
% std_burst_minus = std_burst_timeseries;
% std_mem_minus = std_mem_coef_timeseries;
% 
% addpath('/home/minyoung/data_minyoung/Research/Matlab_function/linspecer')
% C = linspecer(3);
% 
% % burst & memory for minus
% figure2 = figure(2);
% set(gcf,'color','w')
% 
% subplot(1,2,1)
% 
% N = 3;
% C = linspecer(N);
% h1 = plot(mean_burst_minus(:,1),'-','color',C(3,:),'Linewidth',1.5);
% hold on
% h2 = plot(mean_burst_minus(:,2),'-','color',C(1,:),'Linewidth',1.5);
% hold on
% h3 = plot(mean_burst_minus(:,3),'-','color',C(2,:),'Linewidth',1.5);
% hold on
% plot([22 22],[.115 .475],'k--')
% hold on
% ylim([.115 .475])
% 
% hold on
% ciplot(mean_burst_minus(:,1)-cif*std_burst_minus(:,1), mean_burst_minus(:,1)+cif*std_burst_minus(:,1), 1:160, C(3,:))
% hold on
% hold on
% ciplot(mean_burst_minus(:,2)-cif*std_burst_minus(:,2), mean_burst_minus(:,2)+cif*std_burst_minus(:,2), 1:160, C(1,:))
% hold on
% hold on
% ciplot(mean_burst_minus(:,3)-cif*std_burst_minus(:,3), mean_burst_minus(:,3)+cif*std_burst_minus(:,3), 1:160, C(2,:))
% hold on
% alpha(.5)
% 
% date_seq = load('./date/date.txt');
% date_seq = date_seq(win_size:end);
% date_str = cell(size(mean_burst_minus,1),1);
% for i=1:size(mean_burst_minus,1)
%     date_str{i,1} = num2str(date_seq(i,1));
% end
% date_first = [20080815;20080915;20081015;20081114;20081215;20090115;20090216;20090316];
% for i=1:length(date_first)
%     date_first_ind(i,1) = find(date_seq == date_first(i));
% end
% date_str = date_str';
% date_str_tick = date_str(date_first_ind);
% 
% set(gca,'xtick',date_first_ind)
% set(gca,'xticklabel',date_str_tick)
% ax = ancestor(h3,'axes');
% xrule = ax.XAxis;
% xrule.FontSize = 8;
% % ylim([0.14 0.39])
% 
% xtickangle(35)
% xlim([1 160])
% xlabel('Date')
% ylabel('B_1')
% lh = legend([h1 h2 h3],{'Type Zero, Negative','Type A, Negative','Type B, Negative'});
% legend boxoff
% set(lh, 'Position', [-.029 .26 .5 .01])
% % legend('Zero Minus','A Minus','B Minus','Zero Plus','A Plus','B Plus')
% % grid on
% % grid minor
% box on
% text(-0.1,1.1,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
% 
% 
% subplot(1,2,2)
% h1 = plot(mean_mem_minus(:,1),'-','color',C(3,:),'Linewidth',1.5);
% hold on
% h2 = plot(mean_mem_minus(:,2),'-','color',C(1,:),'Linewidth',1.5);
% hold on
% h3 = plot(mean_mem_minus(:,3),'-','color',C(2,:),'Linewidth',1.5);
% hold on
% plot([22 22],[-.045 .26],'k--')
% hold on
% ylim([-.045 .26])
% 
% hold on
% ciplot(mean_mem_minus(:,1)-cif*std_mem_minus(:,1), mean_mem_minus(:,1)+cif*std_mem_minus(:,1), 1:160, C(3,:))
% hold on
% hold on
% ciplot(mean_mem_minus(:,2)-cif*std_mem_minus(:,2), mean_mem_minus(:,2)+cif*std_mem_minus(:,2), 1:160, C(1,:))
% hold on
% hold on
% ciplot(mean_mem_minus(:,3)-cif*std_mem_minus(:,3), mean_mem_minus(:,3)+cif*std_mem_minus(:,3), 1:160, C(2,:))
% hold on
% alpha(.5)
% 
% date_seq = load('./date/date.txt');
% date_seq = date_seq(win_size:end);
% date_str = cell(size(mean_burst_minus,1),1);
% for i=1:size(mean_burst_minus,1)
%     date_str{i,1} = num2str(date_seq(i,1));
% end
% date_first = [20080815;20080915;20081015;20081114;20081215;20090115;20090216;20090316];
% for i=1:length(date_first)
%     date_first_ind(i,1) = find(date_seq == date_first(i));
% end
% date_str = date_str';
% date_str_tick = date_str(date_first_ind);
% 
% set(gca,'xtick',date_first_ind)
% set(gca,'xticklabel',date_str_tick)
% ax = ancestor(h3,'axes');
% xrule = ax.XAxis;
% xrule.FontSize = 8;
% xtickangle(35)
% xlim([1 160])
% xlabel('Date')
% ylabel('M')
% lh = legend([h1 h2 h3],{'Type Zero, Negative','Type A, Negative','Type B, Negative'});
% legend boxoff
% set(lh, 'Position', [.457 .26 .5 .01])
% % legend('Zero Minus','A Minus','B Minus','Zero Plus','A Plus','B Plus')
% % grid on
% % grid minor
% box on
% text(-0.1,1.1,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
% 
% 
% dx=0.06;
% dy=0.06;
% x=(1-4*dx)/1.9;
% y=(1-4*dy)/1.05;
% dxx=0.01;
% dyy=0.01;
% AxesHandle=findobj(figure2,'Type','axes');
% set(AxesHandle(2),'Position',[dx+0.01,3*dy,x,y]);
% set(AxesHandle(1),'Position',[x+2.6*dx,3*dy,x,y]);
% set(gcf,'Position',[100, 100, 900, 400])




