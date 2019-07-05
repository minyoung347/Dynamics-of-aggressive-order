%%
clear;clc;close all;

iter_set = 1; % tick (in order, refer to )
id = 1; % spread (1: all, 2: small, 3: large)
firm_start = 1; % (13)
firm_end = 65;
len_firm = firm_end - firm_start + 1;
day_start = 32; % [2,89], [10,153]
day_end = 61; % (firm:13, day:31) (firm:1-30, day:32-61, error)
t = 1000;


cd '/home/minyoung/data_minyoung/Research/Price_Impact'
filepath='/home/minyoung/data_minyoung/cm_code/test/data/qth_corr_data/result/LSE_data2';
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/');
load('/home/minyoung/data_minyoung/Research/Price_Impact/date/date.txt')
com_num=load('/home/minyoung/data_minyoung/Research/Price_Impact/top_69.csv');


len_day = day_end - day_start + 1 - 1;

% iter_set=1;

aggregate_size = 1;
gamma_crit = 0;

tolerance = 0.00001; % for comparison of two variable of the double type

iter = 0;
NI = 12;


bid_price_minus_ta{1,1}=cell(NI,1);
bid_price_plus_ta{1,1}=cell(NI,1);
ask_price_minus_ta{1,1}=cell(NI,1);
ask_price_plus_ta{1,1}=cell(NI,1);
spread_minus_ta{1,1}=cell(NI,1);
spread_plus_ta{1,1}=cell(NI,1);
% update_time_gap_minus_ta{1,1}=cell(NI,1);
% update_time_gap_plus_ta{1,1}=cell(NI,1);

[bid_price_minus_ta{1,1}{:}]=deal(zeros(2*t+1,1));
[bid_price_plus_ta{1,1}{:}]=deal(zeros(2*t+1,1));
[ask_price_minus_ta{1,1}{:}]=deal(zeros(2*t+1,1));
[ask_price_plus_ta{1,1}{:}]=deal(zeros(2*t+1,1));
[spread_minus_ta{1,1}{:}]=deal(zeros(2*t+1,1));
[spread_plus_ta{1,1}{:}]=deal(zeros(2*t+1,1));
% [update_time_gap_minus_ta{1,1}{:}]=deal(zeros(2*t+1,1));
% [update_time_gap_plus_ta{1,1}{:}]=deal(zeros(2*t+1,1));


bid_vol_minus_ta{1,1}=cell(NI,5);
bid_vol_plus_ta{1,1}=cell(NI,5);
bid_gap_minus_ta{1,1}=cell(NI,5);
bid_gap_plus_ta{1,1}=cell(NI,5);
ask_vol_minus_ta{1,1}=cell(NI,5);
ask_vol_plus_ta{1,1}=cell(NI,5);
ask_gap_minus_ta{1,1}=cell(NI,5);
ask_gap_plus_ta{1,1}=cell(NI,5);
% flow_minus_ta{1,1}=cell(NI,18);
% flow_plus_ta{1,1}=cell(NI,18);
% bid_vol_price_minus_avg_ta{1,1}=cell(NI,4);
% bid_vol_price_plus_avg_ta{1,1}=cell(NI,4);
% ask_vol_price_minus_avg_ta{1,1}=cell(NI,4);
% ask_vol_price_plus_avg_ta{1,1}=cell(NI,4);

[bid_vol_minus_ta{1,1}{:}]=deal(zeros(2*t+1,1));
[bid_vol_plus_ta{1,1}{:}]=deal(zeros(2*t+1,1));
[bid_gap_minus_ta{1,1}{:}]=deal(zeros(2*t+1,1));
[bid_gap_plus_ta{1,1}{:}]=deal(zeros(2*t+1,1));
[ask_vol_minus_ta{1,1}{:}]=deal(zeros(2*t+1,1));
[ask_vol_plus_ta{1,1}{:}]=deal(zeros(2*t+1,1));
[ask_gap_minus_ta{1,1}{:}]=deal(zeros(2*t+1,1));
[ask_gap_plus_ta{1,1}{:}]=deal(zeros(2*t+1,1));
% [flow_minus_ta{1,1}{:}]=deal(zeros(2*t+1,1));
% [flow_plus_ta{1,1}{:}]=deal(zeros(2*t+1,1));
% [bid_vol_price_minus_avg_ta{1,1}{:}]=deal(zeros(2*t+1,1));
% [bid_vol_price_plus_avg_ta{1,1}{:}]=deal(zeros(2*t+1,1));
% [ask_vol_price_minus_avg_ta{1,1}{:}]=deal(zeros(2*t+1,1));
% [ask_vol_price_plus_avg_ta{1,1}{:}]=deal(zeros(2*t+1,1));

sum_minus_ta = zeros(NI,len_firm);
sum_plus_ta = zeros(NI,len_firm);

num_event_minus = zeros(NI,1);
num_event_plus = zeros(NI,1);

num_event_minus_firm = zeros(NI,len_firm);
num_event_plus_firm = zeros(NI,len_firm);

firm_ind = 1;

for firm=firm_start:firm_end
    ini = 1;
    tick_bin_valid_save=[];
    iter = 0;
    day_ind = 1;
    
    sum_minus = zeros(NI,1);
    sum_plus = zeros(NI,1);
    
    for day=day_start:day_end
        [firm,day]
         if exist(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_AB_sign_%d_%d_%d_%d.mat',com_num(firm),day,id,t),'file') && day~=26 ...
                && exist(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_AB_sign_ind_%d_%d_%d_%d.mat',com_num(firm),day,id,t),'file') ...
                && exist(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_AB_sign_var_%d_%d_%d_%d.mat',com_num(firm),day,id,t),'file') ...
                && exist(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_ABa_sign_%d_%d_%d_%d.mat',com_num(firm),day,id,t),'file') ...
                && exist(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_ABa_sign_ind_%d_%d_%d_%d.mat',com_num(firm),day,id,t),'file') ...
                && exist(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result_over/lmc_flow_event_ver6_AB_over_sign_%d_%d_%d_%d.mat',com_num(firm),day,id,t),'file') ...
                && exist(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result_over/lmc_flow_event_ver6_AB_over_sign_ind_%d_%d_%d_%d.mat',com_num(firm),day,id,t),'file')
            
            load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_AB_sign_%d_%d_%d_%d.mat',com_num(firm),day,id,t))
            load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_AB_sign_ind_%d_%d_%d_%d.mat',com_num(firm),day,id,t))
            load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_AB_sign_var_%d_%d_%d_%d.mat',com_num(firm),day,id,t))

            zero_sign = load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_ABa_sign_%d_%d_%d_%d.mat',com_num(firm),day,id,t));
            zero_sign_ind = load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_ABa_sign_ind_%d_%d_%d_%d.mat',com_num(firm),day,id,t));
            over_sign = load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result_over/lmc_flow_event_ver6_AB_over_sign_%d_%d_%d_%d.mat',com_num(firm),day,id,t));
            over_sign_ind = load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result_over/lmc_flow_event_ver6_AB_over_sign_ind_%d_%d_%d_%d.mat',com_num(firm),day,id,t));
            
            % AB + Zero + Over
            ind_minus{10,1} = zero_sign_ind.ind_minus{1,1};
            ind_plus{10,1} = zero_sign_ind.ind_plus{1,1};
            ind_minus{11,1} = over_sign_ind.ind_minus{1,1};
            ind_plus{11,1} = over_sign_ind.ind_plus{1,1};
            ind_minus{12,1} = over_sign_ind.ind_minus{2,1};
            ind_plus{12,1} = over_sign_ind.ind_plus{2,1};
            
            tick_bin_valid = round(100000*tick_bin_valid)/100000; % floating point
            
            for i=1:length(tick_bin_valid)
%                 name_var = 'spread_plus';
%                 eval(sprintf('%s{%d,1}{10,1}',name_var,i)) = [];
%                 eval(sprintf('%s{%d,1}{10,1}',name_var,i)) = eval(sprintf('zero_sign.%s{%d,1}{1,1}',name_var,i));
                spread_plus{i,1}{10,1} = zero_sign.spread_plus{i,1}{1,1};
                spread_plus{i,1}{11,1} = over_sign.spread_plus{i,1}{1,1};
                spread_plus{i,1}{12,1} = over_sign.spread_plus{i,1}{2,1};
                spread_minus{i,1}{10,1} = zero_sign.spread_minus{i,1}{1,1};
                spread_minus{i,1}{11,1} = over_sign.spread_minus{i,1}{1,1};
                spread_minus{i,1}{12,1} = over_sign.spread_minus{i,1}{2,1};
%                 update_time_gap_plus{i,1}{10,1} = zero_sign.update_time_gap_plus{i,1}{1,1};
%                 update_time_gap_plus{i,1}{11,1} = over_sign.update_time_gap_plus{i,1}{1,1};
%                 update_time_gap_plus{i,1}{12,1} = over_sign.update_time_gap_plus{i,1}{2,1};
%                 update_time_gap_minus{i,1}{10,1} = zero_sign.update_time_gap_minus{i,1}{1,1};
%                 update_time_gap_minus{i,1}{11,1} = over_sign.update_time_gap_minus{i,1}{1,1};
%                 update_time_gap_minus{i,1}{12,1} = over_sign.update_time_gap_minus{i,1}{2,1};
                bid_price_plus{i,1}{10,1} = zero_sign.bid_price_plus{i,1}{1,1};
                bid_price_plus{i,1}{11,1} = over_sign.bid_price_plus{i,1}{1,1};
                bid_price_plus{i,1}{12,1} = over_sign.bid_price_plus{i,1}{2,1};
                bid_price_minus{i,1}{10,1} = zero_sign.bid_price_minus{i,1}{1,1};
                bid_price_minus{i,1}{11,1} = over_sign.bid_price_minus{i,1}{1,1};
                bid_price_minus{i,1}{12,1} = over_sign.bid_price_minus{i,1}{2,1};
                ask_price_plus{i,1}{10,1} = zero_sign.ask_price_plus{i,1}{1,1};
                ask_price_plus{i,1}{11,1} = over_sign.ask_price_plus{i,1}{1,1};
                ask_price_plus{i,1}{12,1} = over_sign.ask_price_plus{i,1}{2,1};
                ask_price_minus{i,1}{10,1} = zero_sign.ask_price_minus{i,1}{1,1};
                ask_price_minus{i,1}{11,1} = over_sign.ask_price_minus{i,1}{1,1};
                ask_price_minus{i,1}{12,1} = over_sign.ask_price_minus{i,1}{2,1};
                for ii=1:10
                    bid_vol_plus{i,1}{10,ii} = zero_sign.bid_vol_plus{i,1}{1,ii};
                    bid_vol_plus{i,1}{11,ii} = over_sign.bid_vol_plus{i,1}{1,ii};
                    bid_vol_plus{i,1}{12,ii} = over_sign.bid_vol_plus{i,1}{2,ii};
                    bid_vol_minus{i,1}{10,ii} = zero_sign.bid_vol_minus{i,1}{1,ii};
                    bid_vol_minus{i,1}{11,ii} = over_sign.bid_vol_minus{i,1}{1,ii};
                    bid_vol_minus{i,1}{12,ii} = over_sign.bid_vol_minus{i,1}{2,ii};
                    ask_vol_plus{i,1}{10,ii} = zero_sign.ask_vol_plus{i,1}{1,ii};
                    ask_vol_plus{i,1}{11,ii} = over_sign.ask_vol_plus{i,1}{1,ii};
                    ask_vol_plus{i,1}{12,ii} = over_sign.ask_vol_plus{i,1}{2,ii};
                    ask_vol_minus{i,1}{10,ii} = zero_sign.ask_vol_minus{i,1}{1,ii};
                    ask_vol_minus{i,1}{11,ii} = over_sign.ask_vol_minus{i,1}{1,ii};
                    ask_vol_minus{i,1}{12,ii} = over_sign.ask_vol_minus{i,1}{2,ii};
                    bid_vol_price_plus{i,1}{10,ii} = zero_sign.bid_vol_price_plus{i,1}{1,ii};
                    bid_vol_price_plus{i,1}{11,ii} = over_sign.bid_vol_price_plus{i,1}{1,ii};
                    bid_vol_price_plus{i,1}{12,ii} = over_sign.bid_vol_price_plus{i,1}{2,ii};
                    bid_vol_price_minus{i,1}{10,ii} = zero_sign.bid_vol_price_minus{i,1}{1,ii};
                    bid_vol_price_minus{i,1}{11,ii} = over_sign.bid_vol_price_minus{i,1}{1,ii};
                    bid_vol_price_minus{i,1}{12,ii} = over_sign.bid_vol_price_minus{i,1}{2,ii};
                    ask_vol_price_plus{i,1}{10,ii} = zero_sign.ask_vol_price_plus{i,1}{1,ii};
                    ask_vol_price_plus{i,1}{11,ii} = over_sign.ask_vol_price_plus{i,1}{1,ii};
                    ask_vol_price_plus{i,1}{12,ii} = over_sign.ask_vol_price_plus{i,1}{2,ii};
                    ask_vol_price_minus{i,1}{10,ii} = zero_sign.ask_vol_price_minus{i,1}{1,ii};
                    ask_vol_price_minus{i,1}{11,ii} = over_sign.ask_vol_price_minus{i,1}{1,ii};
                    ask_vol_price_minus{i,1}{12,ii} = over_sign.ask_vol_price_minus{i,1}{2,ii};
                    bid_gap_plus{i,1}{10,ii} = zero_sign.bid_gap_plus{i,1}{1,ii};
                    bid_gap_plus{i,1}{11,ii} = over_sign.bid_gap_plus{i,1}{1,ii};
                    bid_gap_plus{i,1}{12,ii} = over_sign.bid_gap_plus{i,1}{2,ii};
                    bid_gap_minus{i,1}{10,ii} = zero_sign.bid_gap_minus{i,1}{1,ii};
                    bid_gap_minus{i,1}{11,ii} = over_sign.bid_gap_minus{i,1}{1,ii};
                    bid_gap_minus{i,1}{12,ii} = over_sign.bid_gap_minus{i,1}{2,ii};
                    ask_gap_plus{i,1}{10,ii} = zero_sign.ask_gap_plus{i,1}{1,ii};
                    ask_gap_plus{i,1}{11,ii} = over_sign.ask_gap_plus{i,1}{1,ii};
                    ask_gap_plus{i,1}{12,ii} = over_sign.ask_gap_plus{i,1}{2,ii};
                    ask_gap_minus{i,1}{10,ii} = zero_sign.ask_gap_minus{i,1}{1,ii};
                    ask_gap_minus{i,1}{11,ii} = over_sign.ask_gap_minus{i,1}{1,ii};
                    ask_gap_minus{i,1}{12,ii} = over_sign.ask_gap_minus{i,1}{2,ii};
                end
%                 for ii=1:18
%                     flow_plus{i,1}{10,ii} = zero_sign.flow_plus{i,1}{1,ii};
%                     flow_plus{i,1}{11,ii} = over_sign.flow_plus{i,1}{1,ii};
%                     flow_plus{i,1}{12,ii} = over_sign.flow_plus{i,1}{2,ii};
%                     flow_minus{i,1}{10,ii} = zero_sign.flow_minus{i,1}{1,ii};
%                     flow_minus{i,1}{11,ii} = over_sign.flow_minus{i,1}{1,ii};
%                     flow_minus{i,1}{12,ii} = over_sign.flow_minus{i,1}{2,ii};
%                 end
%                 for ii=1:2*t+1
%                     bid_vol_price_p_avg{i,1}{10,ii} = zero_sign.bid_vol_price_p_avg{i,1}{1,ii};
%                     bid_vol_price_p_avg{i,1}{11,ii} = over_sign.bid_vol_price_p_avg{i,1}{1,ii};
%                     bid_vol_price_p_avg{i,1}{12,ii} = over_sign.bid_vol_price_p_avg{i,1}{2,ii};
%                     bid_vol_price_m_avg{i,1}{10,ii} = zero_sign.bid_vol_price_m_avg{i,1}{1,ii};
%                     bid_vol_price_m_avg{i,1}{11,ii} = over_sign.bid_vol_price_m_avg{i,1}{1,ii};
%                     bid_vol_price_m_avg{i,1}{12,ii} = over_sign.bid_vol_price_m_avg{i,1}{2,ii};
%                     ask_vol_price_p_avg{i,1}{10,ii} = zero_sign.ask_vol_price_p_avg{i,1}{1,ii};
%                     ask_vol_price_p_avg{i,1}{11,ii} = over_sign.ask_vol_price_p_avg{i,1}{1,ii};
%                     ask_vol_price_p_avg{i,1}{12,ii} = over_sign.ask_vol_price_p_avg{i,1}{2,ii};
%                     ask_vol_price_m_avg{i,1}{10,ii} = zero_sign.ask_vol_price_m_avg{i,1}{1,ii};
%                     ask_vol_price_m_avg{i,1}{11,ii} = over_sign.ask_vol_price_m_avg{i,1}{1,ii};
%                     ask_vol_price_m_avg{i,1}{12,ii} = over_sign.ask_vol_price_m_avg{i,1}{2,ii};
%                 end
            end
            clear zero_sign zero_sign_ind over_sign over_sign_ind
            
            % remove tick_bin of small size
            sum_num = zeros(length(tick_bin_valid),1);
            for te=1:length(tick_bin_valid)
                for tee=1:NI
                    sum_num(te,1) = sum_num(te,1) + ask_price_minus{te,1}{tee,1}(1,2) + ask_price_plus{te,1}{tee,1}(1,2); 
                    % error ??
                    % sum_num(te,1) = sum_num(te,1) + sum(ask_price_minus{te,1}{tee,1}(:,2));
                end
            end
            
            del_ind = find(sum_num < 100);
            save_ind = find(sum_num >= 100);
            
            if ~isempty(del_ind)
                fid = fopen('to_z15o_jump.txt','a');
                fprintf(fid,'%.0f %.0f %f %.0f \n',firm,day,tick_bin_valid(del_ind),sum_num(del_ind));
                fclose(fid);
            end
            
            tick_bin_valid(del_ind,:)=[];
            ask_gap_minus(del_ind,:)=[];
            ask_gap_plus(del_ind,:)=[];
            ask_price_minus(del_ind,:)=[];
            ask_price_plus(del_ind,:)=[];
            ask_vol_minus(del_ind,:)=[];
            ask_vol_plus(del_ind,:)=[];
            ask_vol_price_m_avg(del_ind,:)=[];
            ask_vol_price_minus(del_ind,:)=[];
            ask_vol_price_p_avg(del_ind,:)=[];
            ask_vol_price_plus(del_ind,:)=[];
            bid_gap_minus(del_ind,:)=[];
            bid_gap_plus(del_ind,:)=[];
            bid_price_minus(del_ind,:)=[];
            bid_price_plus(del_ind,:)=[];
            bid_vol_minus(del_ind,:)=[];
            bid_vol_plus(del_ind,:)=[];
%             bid_vol_price_m_avg(del_ind,:)=[];
            bid_vol_price_minus(del_ind,:)=[];
%             bid_vol_price_p_avg(del_ind,:)=[];
            bid_vol_price_plus(del_ind,:)=[];
%             flow_minus(del_ind,:)=[];
%             flow_plus(del_ind,:)=[];
            spread_minus(del_ind,:)=[];
            spread_plus(del_ind,:)=[];
%             update_time_gap_minus(del_ind,:)=[];
%             update_time_gap_plus(del_ind,:)=[];
            for i=1:NI
                ind_minus{i,1}(del_ind,:)=[];
                ind_plus{i,1}(del_ind,:)=[];
            end
    
%             % change data form
%             for i=1:length(tick_bin_valid)
%                 for ii=1:NI
%                     for iii=1:2*t+1
%                         for iiii=1:4
%                             bid_vol_price_minus_avg{i,1}{ii,iiii}(iii,:) = bid_vol_price_m_avg{i,1}{ii,iii}(iiii,:);
%                             bid_vol_price_plus_avg{i,1}{ii,iiii}(iii,:) = bid_vol_price_p_avg{i,1}{ii,iii}(iiii,:);
%                             ask_vol_price_minus_avg{i,1}{ii,iiii}(iii,:) = ask_vol_price_m_avg{i,1}{ii,iii}(iiii,:);
%                             ask_vol_price_plus_avg{i,1}{ii,iiii}(iii,:) = ask_vol_price_p_avg{i,1}{ii,iii}(iiii,:);
%                             if isnan(bid_vol_price_minus_avg{i,1}{ii,iiii}(iii,1))
%                                 bid_vol_price_minus_avg{i,1}{ii,iiii}(iii,1) = 0;
%                             end
%                             if isnan(bid_vol_price_plus_avg{i,1}{ii,iiii}(iii,1))
%                                 bid_vol_price_plus_avg{i,1}{ii,iiii}(iii,1) = 0;
%                             end
%                             if isnan(ask_vol_price_minus_avg{i,1}{ii,iiii}(iii,1))
%                                 ask_vol_price_minus_avg{i,1}{ii,iiii}(iii,1) = 0;
%                             end
%                             if isnan(ask_vol_price_plus_avg{i,1}{ii,iiii}(iii,1))
%                                 ask_vol_price_plus_avg{i,1}{ii,iiii}(iii,1) = 0;
%                             end
%                         end
%                     end
%                 end
%             end
%             clear bid_vol_price_m_avg bid_vol_price_p_avg ask_vol_price_m_avg ask_vol_price_p_avg
            
            % comparing tick_bin_valid with tick_bin_valid_save (tolerance)
            % ch = 0;
            check_new = length(tick_bin_valid);
            for ii=1:length(tick_bin_valid_save)
                for i=1:length(tick_bin_valid)
                    if find(abs(tick_bin_valid_save(ii) - tick_bin_valid(i)) < tolerance)
%                     if find(tick_bin_vacelllid_save(ii)>tick_bin_valid(i)*(1-tolerance) & tick_bin_valid_save(ii)<tick_bin_valid(i)*(1+tolerance))
                        tick_bin_valid(i)=tick_bin_valid_save(ii);
                        % ch = ch+1;
                        check_new = check_new - 1;
                    end
                end
            end
            
            
            if ini==1 % initiating *_total for each firm
                for iter=1:length(tick_bin_valid)
                    bid_price_minus_total{iter,1}=cell(NI,1);
                    bid_price_plus_total{iter,1}=cell(NI,1);
                    ask_price_minus_total{iter,1}=cell(NI,1);
                    ask_price_plus_total{iter,1}=cell(NI,1);
                    spread_minus_total{iter,1}=cell(NI,1);
                    spread_plus_total{iter,1}=cell(NI,1);
%                     update_time_gap_minus_total{iter,1}=cell(NI,1);
%                     update_time_gap_plus_total{iter,1}=cell(NI,1);

                    [bid_price_minus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [bid_price_plus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [ask_price_minus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [ask_price_plus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [spread_minus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [spread_plus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
%                     [update_time_gap_minus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
%                     [update_time_gap_plus_total{iter,1}{:}]=deal(zeros(2*t+1,1));

                    bid_vol_minus_total{iter,1}=cell(NI,5);
                    bid_vol_plus_total{iter,1}=cell(NI,5);
                    bid_gap_minus_total{iter,1}=cell(NI,5);
                    bid_gap_plus_total{iter,1}=cell(NI,5);
                    ask_vol_minus_total{iter,1}=cell(NI,5);
                    ask_vol_plus_total{iter,1}=cell(NI,5);
                    ask_gap_minus_total{iter,1}=cell(NI,5);
                    ask_gap_plus_total{iter,1}=cell(NI,5);
%                     flow_minus_total{iter,1}=cell(NI,18);
%                     flow_plus_total{iter,1}=cell(NI,18);
%                     bid_vol_price_minus_avg_total{iter,1}=cell(NI,4);
%                     bid_vol_price_plus_avg_total{iter,1}=cell(NI,4);
%                     ask_vol_price_minus_avg_total{iter,1}=cell(NI,4);
%                     ask_vol_price_plus_avg_total{iter,1}=cell(NI,4);

                    [bid_vol_minus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [bid_vol_plus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [bid_gap_minus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [bid_gap_plus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [ask_vol_minus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [ask_vol_plus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [ask_gap_minus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [ask_gap_plus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
%                     [flow_minus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
%                     [flow_plus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
%                     [bid_vol_price_minus_avg_total{iter,1}{:}]=deal(zeros(2*t+1,1));
%                     [bid_vol_price_plus_avg_total{iter,1}{:}]=deal(zeros(2*t+1,1));
%                     [ask_vol_price_minus_avg_total{iter,1}{:}]=deal(zeros(2*t+1,1));
%                     [ask_vol_price_plus_avg_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    
                    bid_price_minus_num{iter,1}=zeros(NI,1);
                    bid_price_plus_num{iter,1}=zeros(NI,1);
                    ask_price_minus_num{iter,1}=zeros(NI,1);
                    ask_price_plus_num{iter,1}=zeros(NI,1);
                    spread_minus_num{iter,1}=zeros(NI,1);
                    spread_plus_num{iter,1}=zeros(NI,1);
%                     update_time_gap_minus_num{iter,1}=zeros(NI,1);
%                     update_time_gap_plus_num{iter,1}=zeros(NI,1);

                    bid_vol_minus_num{iter,1}=zeros(NI,5);
                    bid_vol_plus_num{iter,1}=zeros(NI,5);
                    bid_gap_minus_num{iter,1}=zeros(NI,5);
                    bid_gap_plus_num{iter,1}=zeros(NI,5);
                    ask_vol_minus_num{iter,1}=zeros(NI,5);
                    ask_vol_plus_num{iter,1}=zeros(NI,5);
                    ask_gap_minus_num{iter,1}=zeros(NI,5);
                    ask_gap_plus_num{iter,1}=zeros(NI,5);
%                     flow_minus_num{iter,1}=zeros(NI,18);
%                     flow_plus_num{iter,1}=zeros(NI,18);
%                     bid_vol_price_minus_avg_num{iter,1}=zeros(NI,4);
%                     bid_vol_price_plus_avg_num{iter,1}=zeros(NI,4);
%                     ask_vol_price_minus_avg_num{iter,1}=zeros(NI,4);
%                     ask_vol_price_plus_avg_num{iter,1}=zeros(NI,4);
                    
                    %
                    bid_price_minus_num{iter,2}=zeros(NI,1);
                    bid_price_plus_num{iter,2}=zeros(NI,1);
                    ask_price_minus_num{iter,2}=zeros(NI,1);
                    ask_price_plus_num{iter,2}=zeros(NI,1);
                    spread_minus_num{iter,2}=zeros(NI,1);
                    spread_plus_num{iter,2}=zeros(NI,1);
%                     update_time_gap_minus_num{iter,2}=zeros(NI,1);
%                     update_time_gap_plus_num{iter,2}=zeros(NI,1);

                    bid_vol_minus_num{iter,2}=zeros(NI,5);
                    bid_vol_plus_num{iter,2}=zeros(NI,5);
                    bid_gap_minus_num{iter,2}=zeros(NI,5);
                    bid_gap_plus_num{iter,2}=zeros(NI,5);
                    ask_vol_minus_num{iter,2}=zeros(NI,5);
                    ask_vol_plus_num{iter,2}=zeros(NI,5);
                    ask_gap_minus_num{iter,2}=zeros(NI,5);
                    ask_gap_plus_num{iter,2}=zeros(NI,5);
%                     flow_minus_num{iter,2}=zeros(NI,18);
%                     flow_plus_num{iter,2}=zeros(NI,18);
%                     bid_vol_price_minus_avg_num{iter,2}=zeros(NI,4);
%                     bid_vol_price_plus_avg_num{iter,2}=zeros(NI,4);
%                     ask_vol_price_minus_avg_num{iter,2}=zeros(NI,4);
%                     ask_vol_price_plus_avg_num{iter,2}=zeros(NI,4);
                end
                ini=0;
            else % after initiating, 'new' tick_bin_valid comes
%                 if ch~=length(tick_bin_valid)
                for cn = 1:check_new
                    iter=iter+1;
                    bid_price_minus_total{iter,1}=cell(NI,1);
                    bid_price_plus_total{iter,1}=cell(NI,1);
                    ask_price_minus_total{iter,1}=cell(NI,1);
                    ask_price_plus_total{iter,1}=cell(NI,1);
                    spread_minus_total{iter,1}=cell(NI,1);
                    spread_plus_total{iter,1}=cell(NI,1);
%                     update_time_gap_minus_total{iter,1}=cell(NI,1);
%                     update_time_gap_plus_total{iter,1}=cell(NI,1);

                    [bid_price_minus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [bid_price_plus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [ask_price_minus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [ask_price_plus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [spread_minus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [spread_plus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
%                     [update_time_gap_minus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
%                     [update_time_gap_plus_total{iter,1}{:}]=deal(zeros(2*t+1,1));

                    bid_vol_minus_total{iter,1}=cell(NI,5);
                    bid_vol_plus_total{iter,1}=cell(NI,5);
                    bid_gap_minus_total{iter,1}=cell(NI,5);
                    bid_gap_plus_total{iter,1}=cell(NI,5);
                    ask_vol_minus_total{iter,1}=cell(NI,5);
                    ask_vol_plus_total{iter,1}=cell(NI,5);
                    ask_gap_minus_total{iter,1}=cell(NI,5);
                    ask_gap_plus_total{iter,1}=cell(NI,5);
%                     flow_minus_total{iter,1}=cell(NI,18);
%                     flow_plus_total{iter,1}=cell(NI,18);
%                     bid_vol_price_minus_avg_total{iter,1}=cell(NI,4);
%                     bid_vol_price_plus_avg_total{iter,1}=cell(NI,4);
%                     ask_vol_price_minus_avg_total{iter,1}=cell(NI,4);
%                     ask_vol_price_plus_avg_total{iter,1}=cell(NI,4);

                    [bid_vol_minus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [bid_vol_plus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [bid_gap_minus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [bid_gap_plus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [ask_vol_minus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [ask_vol_plus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [ask_gap_minus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    [ask_gap_plus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
%                     [flow_minus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
%                     [flow_plus_total{iter,1}{:}]=deal(zeros(2*t+1,1));
%                     [bid_vol_price_minus_avg_total{iter,1}{:}]=deal(zeros(2*t+1,1));
%                     [bid_vol_price_plus_avg_total{iter,1}{:}]=deal(zeros(2*t+1,1));
%                     [ask_vol_price_minus_avg_total{iter,1}{:}]=deal(zeros(2*t+1,1));
%                     [ask_vol_price_plus_avg_total{iter,1}{:}]=deal(zeros(2*t+1,1));
                    
                    bid_price_minus_num{iter,1}=zeros(NI,1);
                    bid_price_plus_num{iter,1}=zeros(NI,1);
                    ask_price_minus_num{iter,1}=zeros(NI,1);
                    ask_price_plus_num{iter,1}=zeros(NI,1);
                    spread_minus_num{iter,1}=zeros(NI,1);
                    spread_plus_num{iter,1}=zeros(NI,1);
%                     update_time_gap_minus_num{iter,1}=zeros(NI,1);
%                     update_time_gap_plus_num{iter,1}=zeros(NI,1);

                    bid_vol_minus_num{iter,1}=zeros(NI,5);
                    bid_vol_plus_num{iter,1}=zeros(NI,5);
                    bid_gap_minus_num{iter,1}=zeros(NI,5);
                    bid_gap_plus_num{iter,1}=zeros(NI,5);
                    ask_vol_minus_num{iter,1}=zeros(NI,5);
                    ask_vol_plus_num{iter,1}=zeros(NI,5);
                    ask_gap_minus_num{iter,1}=zeros(NI,5);
                    ask_gap_plus_num{iter,1}=zeros(NI,5);
%                     flow_minus_num{iter,1}=zeros(NI,18);
%                     flow_plus_num{iter,1}=zeros(NI,18);
%                     bid_vol_price_minus_avg_num{iter,1}=zeros(NI,4);
%                     bid_vol_price_plus_avg_num{iter,1}=zeros(NI,4);
%                     ask_vol_price_minus_avg_num{iter,1}=zeros(NI,4);
%                     ask_vol_price_plus_avg_num{iter,1}=zeros(NI,4);
                    
                    %
                    bid_price_minus_num{iter,2}=zeros(NI,1);
                    bid_price_plus_num{iter,2}=zeros(NI,1);
                    ask_price_minus_num{iter,2}=zeros(NI,1);
                    ask_price_plus_num{iter,2}=zeros(NI,1);
                    spread_minus_num{iter,2}=zeros(NI,1);
                    spread_plus_num{iter,2}=zeros(NI,1);
%                     update_time_gap_minus_num{iter,2}=zeros(NI,1);
%                     update_time_gap_plus_num{iter,2}=zeros(NI,1);

                    bid_vol_minus_num{iter,2}=zeros(NI,5);
                    bid_vol_plus_num{iter,2}=zeros(NI,5);
                    bid_gap_minus_num{iter,2}=zeros(NI,5);
                    bid_gap_plus_num{iter,2}=zeros(NI,5);
                    ask_vol_minus_num{iter,2}=zeros(NI,5);
                    ask_vol_plus_num{iter,2}=zeros(NI,5);
                    ask_gap_minus_num{iter,2}=zeros(NI,5);
                    ask_gap_plus_num{iter,2}=zeros(NI,5);
%                     flow_minus_num{iter,2}=zeros(NI,18);
%                     flow_plus_num{iter,2}=zeros(NI,18);
%                     bid_vol_price_minus_avg_num{iter,2}=zeros(NI,4);
%                     bid_vol_price_plus_avg_num{iter,2}=zeros(NI,4);
%                     ask_vol_price_minus_avg_num{iter,2}=zeros(NI,4);
%                     ask_vol_price_plus_avg_num{iter,2}=zeros(NI,4);
                end
            end

            % expand tick_bin_valid_save
            tick_bin_valid_save = [tick_bin_valid_save;tick_bin_valid];
            [dummy ii] = unique(tick_bin_valid_save,'first');
            tick_bin_valid_save = tick_bin_valid_save(sort(ii));

            % accumulate inter-day data to _total (firm level)
            for ti=1:length(tick_bin_valid)
                ind=find(tick_bin_valid_save==tick_bin_valid(ti,1));
                for i=1:NI
                    for ii=1:5
                        bid_vol_minus_total{ind,1}{i,ii} = bid_vol_minus_total{ind,1}{i,ii} + bid_vol_minus{ti,1}{i,ii}(:,1).*bid_vol_minus{ti,1}{i,ii}(:,2);
                        bid_vol_plus_total{ind,1}{i,ii} = bid_vol_plus_total{ind,1}{i,ii} + bid_vol_plus{ti,1}{i,ii}(:,1).*bid_vol_plus{ti,1}{i,ii}(:,2);
                        bid_gap_minus_total{ind,1}{i,ii} = bid_gap_minus_total{ind,1}{i,ii} + bid_gap_minus{ti,1}{i,ii}(:,1).*bid_gap_minus{ti,1}{i,ii}(:,2);
                        bid_gap_plus_total{ind,1}{i,ii} = bid_gap_plus_total{ind,1}{i,ii} + bid_gap_plus{ti,1}{i,ii}(:,1).*bid_gap_plus{ti,1}{i,ii}(:,2);
                        ask_vol_minus_total{ind,1}{i,ii} = ask_vol_minus_total{ind,1}{i,ii} + ask_vol_minus{ti,1}{i,ii}(:,1).*ask_vol_minus{ti,1}{i,ii}(:,2);
                        ask_vol_plus_total{ind,1}{i,ii} = ask_vol_plus_total{ind,1}{i,ii} + ask_vol_plus{ti,1}{i,ii}(:,1).*ask_vol_plus{ti,1}{i,ii}(:,2);
                        ask_gap_minus_total{ind,1}{i,ii} = ask_gap_minus_total{ind,1}{i,ii} + ask_gap_minus{ti,1}{i,ii}(:,1).*ask_gap_minus{ti,1}{i,ii}(:,2);
                        ask_gap_plus_total{ind,1}{i,ii} = ask_gap_plus_total{ind,1}{i,ii} + ask_gap_plus{ti,1}{i,ii}(:,1).*ask_gap_plus{ti,1}{i,ii}(:,2);

                        bid_vol_minus_num{ind,1}(i,ii) = bid_vol_minus_num{ind,1}(i,ii) + bid_vol_minus{ti,1}{i,ii}(1,2);
                        bid_vol_plus_num{ind,1}(i,ii) = bid_vol_plus_num{ind,1}(i,ii) + bid_vol_plus{ti,1}{i,ii}(1,2);
                        bid_gap_minus_num{ind,1}(i,ii) = bid_gap_minus_num{ind,1}(i,ii) + bid_gap_minus{ti,1}{i,ii}(1,2);
                        bid_gap_plus_num{ind,1}(i,ii) = bid_gap_plus_num{ind,1}(i,ii) + bid_gap_plus{ti,1}{i,ii}(1,2);
                        ask_vol_minus_num{ind,1}(i,ii) = ask_vol_minus_num{ind,1}(i,ii) + ask_vol_minus{ti,1}{i,ii}(1,2);
                        ask_vol_plus_num{ind,1}(i,ii) = ask_vol_plus_num{ind,1}(i,ii) + ask_vol_plus{ti,1}{i,ii}(1,2);
                        ask_gap_minus_num{ind,1}(i,ii) = ask_gap_minus_num{ind,1}(i,ii) + ask_gap_minus{ti,1}{i,ii}(1,2);
                        ask_gap_plus_num{ind,1}(i,ii) = ask_gap_plus_num{ind,1}(i,ii) + ask_gap_plus{ti,1}{i,ii}(1,2);
                        
                        if bid_vol_minus_num{ind,1}(i,ii) > 0
                            bid_vol_minus_num{ind,2}(i,ii) = bid_vol_minus_num{ind,2}(i,ii) + 1;
                            bid_gap_minus_num{ind,2}(i,ii) = bid_gap_minus_num{ind,2}(i,ii) + 1;
                            ask_vol_minus_num{ind,2}(i,ii) = ask_vol_minus_num{ind,2}(i,ii) + 1;
                            ask_gap_minus_num{ind,2}(i,ii) = ask_gap_minus_num{ind,2}(i,ii) + 1;
                        end
                        if bid_vol_plus_num{ind,2}(i,ii) > 0
                            bid_vol_plus_num{ind,2}(i,ii) = bid_vol_plus_num{ind,2}(i,ii) + 1;
                            bid_gap_plus_num{ind,2}(i,ii) = bid_gap_plus_num{ind,2}(i,ii) + 1;
                            ask_vol_plus_num{ind,2}(i,ii) = ask_vol_plus_num{ind,2}(i,ii) + 1;
                            ask_gap_plus_num{ind,2}(i,ii) = ask_gap_plus_num{ind,2}(i,ii) + 1;
                        end
                    end
%                     for ii=1:4
%                         bid_vol_price_minus_avg_total{ind,1}{i,ii} = bid_vol_price_minus_avg_total{ind,1}{i,ii} + bid_vol_price_minus_avg{ti,1}{i,ii}(:,1).*bid_vol_price_minus_avg{ti,1}{i,ii}(:,2);
%                         bid_vol_price_plus_avg_total{ind,1}{i,ii} = bid_vol_price_plus_avg_total{ind,1}{i,ii} + bid_vol_price_plus_avg{ti,1}{i,ii}(:,1).*bid_vol_price_plus_avg{ti,1}{i,ii}(:,2);
%                         ask_vol_price_minus_avg_total{ind,1}{i,ii} = ask_vol_price_minus_avg_total{ind,1}{i,ii} + ask_vol_price_minus_avg{ti,1}{i,ii}(:,1).*ask_vol_price_minus_avg{ti,1}{i,ii}(:,2);
%                         ask_vol_price_plus_avg_total{ind,1}{i,ii} = ask_vol_price_plus_avg_total{ind,1}{i,ii} + ask_vol_price_plus_avg{ti,1}{i,ii}(:,1).*ask_vol_price_plus_avg{ti,1}{i,ii}(:,2);
%                         
%                         bid_vol_price_minus_avg_num{ind,1}(i,ii) = bid_vol_price_minus_avg_num{ind,1}(i,ii) + bid_vol_price_minus_avg{ti,1}{i,ii}(1,2);
%                         bid_vol_price_plus_avg_num{ind,1}(i,ii) = bid_vol_price_plus_avg_num{ind,1}(i,ii) + bid_vol_price_plus_avg{ti,1}{i,ii}(1,2);
%                         ask_vol_price_minus_avg_num{ind,1}(i,ii) = ask_vol_price_minus_avg_num{ind,1}(i,ii) + ask_vol_price_minus_avg{ti,1}{i,ii}(1,2);
%                         ask_vol_price_plus_avg_num{ind,1}(i,ii) = ask_vol_price_plus_avg_num{ind,1}(i,ii) + ask_vol_price_plus_avg{ti,1}{i,ii}(1,2);
%                         
%                         if bid_vol_price_minus_avg_num{ind,1}(i,ii) > 0
%                             bid_vol_price_minus_avg_num{ind,2}(i,ii) = bid_vol_price_minus_avg_num{ind,2}(i,ii) + 1;
%                             ask_vol_price_minus_avg_num{ind,2}(i,ii) = ask_vol_price_minus_avg_num{ind,2}(i,ii) + 1;
%                         end
%                         if bid_vol_price_plus_avg_num{ind,1}(i,ii) > 0
%                             bid_vol_price_plus_avg_num{ind,2}(i,ii) = bid_vol_price_plus_avg_num{ind,2}(i,ii) + 1;
%                             ask_vol_price_plus_avg_num{ind,2}(i,ii) = ask_vol_price_plus_avg_num{ind,2}(i,ii) + 1;
%                         end
%                     end
%                     for ii=1:18
%                         flow_minus_total{ind,1}{i,ii} = flow_minus_total{ind,1}{i,ii} + flow_minus{ti,1}{i,ii}(:,1).*flow_minus{ti,1}{i,ii}(:,2);
%                         flow_plus_total{ind,1}{i,ii} = flow_plus_total{ind,1}{i,ii} + flow_plus{ti,1}{i,ii}(:,1).*flow_plus{ti,1}{i,ii}(:,2);
%                         
%                         flow_minus_num{ind,1}(i,ii) = flow_minus_num{ind,1}(i,ii) + flow_minus{ti,1}{i,ii}(1,2);
%                         flow_plus_num{ind,1}(i,ii) = flow_plus_num{ind,1}(i,ii) + flow_plus{ti,1}{i,ii}(1,2);
%                         
%                         if flow_minus_num{ind,1}(i,ii) > 0
%                             flow_minus_num{ind,2}(i,ii) = flow_minus_num{ind,2}(i,ii) + 1;
%                         end
%                         if flow_plus_num{ind,1}(i,ii) > 0
%                             flow_plus_num{ind,2}(i,ii) = flow_plus_num{ind,2}(i,ii) + 1;
%                         end
%                     end
                    bid_price_minus_total{ind,1}{i,1} = bid_price_minus_total{ind,1}{i,1} + bid_price_minus{ti,1}{i,1}(:,1).*bid_price_minus{ti,1}{i,1}(:,2);
                    bid_price_plus_total{ind,1}{i,1} = bid_price_plus_total{ind,1}{i,1} + bid_price_plus{ti,1}{i,1}(:,1).*bid_price_plus{ti,1}{i,1}(:,2);
                    ask_price_minus_total{ind,1}{i,1} = ask_price_minus_total{ind,1}{i,1} + ask_price_minus{ti,1}{i,1}(:,1).*ask_price_minus{ti,1}{i,1}(:,2);
                    ask_price_plus_total{ind,1}{i,1} = ask_price_plus_total{ind,1}{i,1} + ask_price_plus{ti,1}{i,1}(:,1).*ask_price_plus{ti,1}{i,1}(:,2);
                    spread_minus_total{ind,1}{i,1} = spread_minus_total{ind,1}{i,1} + spread_minus{ti,1}{i,1}(:,1).*spread_minus{ti,1}{i,1}(:,2);
                    spread_plus_total{ind,1}{i,1} = spread_plus_total{ind,1}{i,1} + spread_plus{ti,1}{i,1}(:,1).*spread_plus{ti,1}{i,1}(:,2);
%                     update_time_gap_minus_total{ind,1}{i,1} = update_time_gap_minus_total{ind,1}{i,1} + update_time_gap_minus{ti,1}{i,1}(:,1).*update_time_gap_minus{ti,1}{i,1}(:,2);
%                     update_time_gap_plus_total{ind,1}{i,1} = update_time_gap_plus_total{ind,1}{i,1} + update_time_gap_plus{ti,1}{i,1}(:,1).*update_time_gap_plus{ti,1}{i,1}(:,2);

                    bid_price_minus_num{ind,1}(i,1) = bid_price_minus_num{ind,1}(i,1) + bid_price_minus{ti,1}{i,1}(1,2);
                    bid_price_plus_num{ind,1}(i,1) = bid_price_plus_num{ind,1}(i,1) + bid_price_plus{ti,1}{i,1}(1,2);
                    ask_price_minus_num{ind,1}(i,1) = ask_price_minus_num{ind,1}(i,1) + ask_price_minus{ti,1}{i,1}(1,2);
                    ask_price_plus_num{ind,1}(i,1) = ask_price_plus_num{ind,1}(i,1) + ask_price_plus{ti,1}{i,1}(1,2);
                    spread_minus_num{ind,1}(i,1) = spread_minus_num{ind,1}(i,1) + spread_minus{ti,1}{i,1}(1,2);
                    spread_plus_num{ind,1}(i,1) = spread_plus_num{ind,1}(i,1) + spread_plus{ti,1}{i,1}(1,2);
%                     update_time_gap_minus_num{ind,1}(i,1) = update_time_gap_minus_num{ind,1}(i,1) + update_time_gap_minus{ti,1}{i,1}(1,2);
%                     update_time_gap_plus_num{ind,1}(i,1) = update_time_gap_plus_num{ind,1}(i,1) + update_time_gap_plus{ti,1}{i,1}(1,2);
                    
                    if bid_price_minus_num{ind,1}(i,1) > 0
                        bid_price_minus_num{ind,2}(i,1) = bid_price_minus_num{ind,2}(i,1) + 1;
                        ask_price_minus_num{ind,2}(i,1) = ask_price_minus_num{ind,2}(i,1) + 1;
                        spread_minus_num{ind,2}(i,1) = spread_minus_num{ind,2}(i,1) + 1;
%                         update_time_gap_minus_num{ind,2}(i,1) = update_time_gap_minus_num{ind,2}(i,1) + 1;
                    end
                    if bid_price_plus_num{ind,1}(i,1) > 0
                        bid_price_plus_num{ind,2}(i,1) = bid_price_plus_num{ind,2}(i,1) + 1;
                        ask_price_plus_num{ind,2}(i,1) = ask_price_plus_num{ind,2}(i,1) + 1;
                        spread_plus_num{ind,2}(i,1) = spread_plus_num{ind,2}(i,1) + 1;
%                         update_time_gap_plus_num{ind,2}(i,1) = update_time_gap_plus_num{ind,2}(i,1) + 1;
                    end
                end
            end



            clear ind_minus ind_plus
    %         clear ind_minus_0 ind_minus_1 ind_minus_2 ind_minus_3 ind_minus_4 ind_minus_5 ind_minus_6 ind_plus_0 ind_plus_1 ind_plus_2 ind_plus_3 ind_plus_4 ind_plus_5 ind_plus_6
    %         clear ind_minus_7 ind_minus_8 ind_minus_9 ind_minus_10 ind_minus_11 ind_plus_7 ind_plus_8 ind_plus_9 ind_plus_10 ind_plus_11
    %         clear ind_minus_12 ind_minus_13 ind_minus_14 ind_minus_15 ind_minus_16 ind_minus_17 ind_minus_18 ind_minus_19
    %         clear ind_plus_12 ind_plus_13 ind_plus_14 ind_plus_15 ind_plus_16 ind_plus_17 ind_plus_18 ind_plus_19
            clear ask_gap_minus ask_gap_plus ask_price_minus ask_price_plus ask_vol_minus ask_vol_plus
            clear bid_gap_minus bid_gap_plus bid_price_minus bid_price_plus bid_vol_minus bid_vol_plus spread_minus spread_plus flow_minus flow_plus
            fid = fopen('to_z15o_pass.txt','a');
            fprintf(fid,'%.0f %.0f\n',firm,day);
            fclose(fid);
        else
            fid = fopen('to_z15o_fail.txt','a');
            fprintf(fid,'%.0f %.0f\n',firm,day);
            fclose(fid);
        end
        clear bid_vol_price_minus_avg bid_vol_price_plus_avg ask_vol_price_minus_avg ask_vol_price_plus_avg
        clear aa ask ask_di ask_gap ask_gap_event ask_price ask_price_event ask_return ask_vol ask_vol_event
        clear bid bid_di bid_gap bid_gap_event bid_price bid_price_event bid_return bid_vol bid_vol_event
        clear cancel flow lc_minus limit market orderflow price price_return spread spread_a spread_b t_time
        clear time_event time_gap flow trade_id update_time_gap_minus update_time_gap_plus
        clear ask_gap_event_ind ask_price_event_ind ask_return_ind ask_vol_event_ind bid_gap_event_ind bid_price_event_ind
        clear bid_return_ind bid_vol_event_ind cancel_ind lc_minus_ind limit_ind market_ind orderflow_ind
        clear price_ind price_return_ind spread_ind t_time_ind time_event_ind time_gap_ind flow_ind
% day
    day_ind = day_ind + 1;
    end

%     for iter=1:size(tick_bin_valid_save,1)
%         spread_plus_total{iter,1}{i,1} = spread_plus_total{iter,1}{i,1}/spread_plus_num{iter,1}(i,1);
%     
%     end
    % /num
    
    % divided by number of event and tick_bin_valid
    for ti=1:length(tick_bin_valid_save)
        for i=1:NI
            if bid_price_minus_num{ti,1}(i,1) > 0
                for ii=1:5
                    bid_vol_minus_total{ti,1}{i,ii} = bid_vol_minus_total{ti,1}{i,ii};
                    bid_gap_minus_total{ti,1}{i,ii} = bid_gap_minus_total{ti,1}{i,ii} / tick_bin_valid_save(ti);
                    ask_vol_minus_total{ti,1}{i,ii} = ask_vol_minus_total{ti,1}{i,ii};
                    ask_gap_minus_total{ti,1}{i,ii} = ask_gap_minus_total{ti,1}{i,ii} / tick_bin_valid_save(ti);
                end
%                 for ii=1:4
%                     bid_vol_price_minus_avg_total{ti,1}{i,ii} = bid_vol_price_minus_avg_total{ti,1}{i,ii} / tick_bin_valid_save(ti);
%                     ask_vol_price_minus_avg_total{ti,1}{i,ii} = ask_vol_price_minus_avg_total{ti,1}{i,ii} / tick_bin_valid_save(ti);
%                 end
%                 for ii=1:18
%                     flow_minus_total{ti,1}{i,ii} = flow_minus_total{ti,1}{i,ii};
%                 end
                bid_price_minus_total{ti,1}{i,1} = bid_price_minus_total{ti,1}{i,1} / tick_bin_valid_save(ti);
                ask_price_minus_total{ti,1}{i,1} = ask_price_minus_total{ti,1}{i,1} / tick_bin_valid_save(ti);
                spread_minus_total{ti,1}{i,1} = spread_minus_total{ti,1}{i,1}; % spread is considerd in [calculation_event_ver6.m]
%                 update_time_gap_minus_total{ti,1}{i,1} = update_time_gap_minus_total{ti,1}{i,1};
            end
            
            if bid_price_plus_num{ti,1}(i,1) > 0
                for ii=1:5
                    bid_vol_plus_total{ti,1}{i,ii} = bid_vol_plus_total{ti,1}{i,ii};
                    bid_gap_plus_total{ti,1}{i,ii} = bid_gap_plus_total{ti,1}{i,ii} / tick_bin_valid_save(ti);
                    ask_vol_plus_total{ti,1}{i,ii} = ask_vol_plus_total{ti,1}{i,ii};
                    ask_gap_plus_total{ti,1}{i,ii} = ask_gap_plus_total{ti,1}{i,ii} / tick_bin_valid_save(ti);
                end
%                 for ii=1:4
%                     bid_vol_price_plus_avg_total{ti,1}{i,ii} = bid_vol_price_plus_avg_total{ti,1}{i,ii} / tick_bin_valid_save(ti);
%                     ask_vol_price_plus_avg_total{ti,1}{i,ii} = ask_vol_price_plus_avg_total{ti,1}{i,ii} / tick_bin_valid_save(ti);
%                 end
%                 for ii=1:18
%                     flow_plus_total{ti,1}{i,ii} = flow_plus_total{ti,1}{i,ii};
%                 end
                bid_price_plus_total{ti,1}{i,1} = bid_price_plus_total{ti,1}{i,1} / tick_bin_valid_save(ti);
                ask_price_plus_total{ti,1}{i,1} = ask_price_plus_total{ti,1}{i,1} / tick_bin_valid_save(ti);
                spread_plus_total{ti,1}{i,1} = spread_plus_total{ti,1}{i,1};
%                 update_time_gap_plus_total{ti,1}{i,1} = update_time_gap_plus_total{ti,1}{i,1};
            end
        end
    end
    
    for ti=1:length(tick_bin_valid_save)
        for i=1:NI
            num_event_minus(i,1) = num_event_minus(i,1) + bid_price_minus_num{ti,1}(i,1);
            num_event_plus(i,1) = num_event_plus(i,1) + bid_price_plus_num{ti,1}(i,1);
            
            num_event_minus_firm(i,firm_ind) = num_event_minus_firm(i,firm_ind) + bid_price_minus_num{ti,1}(i,1);
            num_event_plus_firm(i,firm_ind) = num_event_plus_firm(i,firm_ind) + bid_price_plus_num{ti,1}(i,1);
        end
    end
    
%     % normalize
%     for ti=1:length(tick_bin_valid_save)
%         for i=1:NI
%             for ii=1:5
% %                 a = bid_vol_minus_total{ti,1}{i,ii}; bid_vol_minus_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
% %                 a = bid_vol_plus_total{ti,1}{i,ii}; bid_vol_plus_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
% %                 a = bid_vol_plus_total{ti,1}{i,ii}; bid_vol_plus_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
%                 a = bid_gap_minus_total{ti,1}{i,ii}; bid_gap_minus_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
%                 a = bid_gap_plus_total{ti,1}{i,ii}; bid_gap_plus_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
% %                 a = ask_vol_minus_total{ti,1}{i,ii}; ask_vol_minus_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
% %                 a = ask_vol_plus_total{ti,1}{i,ii}; ask_vol_plus_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
%                 a = ask_gap_minus_total{ti,1}{i,ii}; ask_gap_minus_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
%                 a = ask_gap_plus_total{ti,1}{i,ii}; ask_gap_plus_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
%                 
%                 
% %                 a = bid_vol_minus_total{ti,1}{i,ii}; bid_vol_minus_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
% %                 a = bid_vol_plus_total{ti,1}{i,ii}; bid_vol_plus_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
% % %                 a = bid_vol_plus_total{ti,1}{i,ii}; bid_vol_plus_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
% %                 a = bid_gap_minus_total{ti,1}{i,ii}; bid_gap_minus_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
% %                 a = bid_gap_plus_total{ti,1}{i,ii}; bid_gap_plus_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
% %                 a = ask_vol_minus_total{ti,1}{i,ii}; ask_vol_minus_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
% %                 a = ask_vol_plus_total{ti,1}{i,ii}; ask_vol_plus_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
% %                 a = ask_gap_minus_total{ti,1}{i,ii}; ask_gap_minus_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
% %                 a = ask_gap_plus_total{ti,1}{i,ii}; ask_gap_plus_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
%             end
%             for ii=1:4
% %                 a = bid_vol_price_minus_avg_total{ti,1}{i,ii}; bid_vol_price_minus_avg_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
% %                 a = bid_vol_price_plus_avg_total{ti,1}{i,ii}; bid_vol_price_plus_avg_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
% %                 a = ask_vol_price_minus_avg_total{ti,1}{i,ii}; ask_vol_price_minus_avg_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
% %                 a = ask_vol_price_plus_avg_total{ti,1}{i,ii}; ask_vol_price_plus_avg_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
%                 a = bid_vol_price_minus_avg_total{ti,1}{i,ii}; bid_vol_price_minus_avg_total{ti,1}{i,ii} = a;
%                 a = bid_vol_price_plus_avg_total{ti,1}{i,ii}; bid_vol_price_plus_avg_total{ti,1}{i,ii} = a;
%                 a = ask_vol_price_minus_avg_total{ti,1}{i,ii}; ask_vol_price_minus_avg_total{ti,1}{i,ii} = a;
%                 a = ask_vol_price_plus_avg_total{ti,1}{i,ii}; ask_vol_price_plus_avg_total{ti,1}{i,ii} = a;
%             end
%             for ii=1:18
%                 a = flow_minus_total{ti,1}{i,ii}; flow_minus_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
%                 a = flow_plus_total{ti,1}{i,ii}; flow_plus_total{ti,1}{i,ii} = (a-a(t+1))/1 + 100;
%             end
%             
% % % %             refer_use =  bid_price_minus_total{iter,1}{i,1}(t+1);
% % % %             b = [];
% % % %             for ai = 1:NI
% % % %                 b = [b; bid_price_minus_total{ti,1}{ai,1}];
% % % %             end
% % % %             a = bid_price_minus_total{ti,1}{i,1}; bid_price_minus_total{ti,1}{i,1} = (a-refer_use);
% % % %             
% % % %             b = [];
% % % %             for ai = 1:NI
% % % %                 b = [b; ask_price_minus_total{ti,1}{ai,1}];
% % % %             end
% % % %             a = ask_price_minus_total{ti,1}{i,1}; ask_price_minus_total{ti,1}{i,1} = (a-refer_use);
% % % %             
% % % %             
% % % %             
% % % %             refer_use =  bid_price_plus_total{iter,1}{i,1}(t+1);
% % % %             b = [];
% % % %             for ai = 1:NI
% % % %                 b = [b; bid_price_plus_total{ti,1}{ai,1}];
% % % %             end
% % % %             a = bid_price_plus_total{ti,1}{i,1}; bid_price_plus_total{ti,1}{i,1} = (a-refer_use);
% % % %             
% % % %             b = [];
% % % %             for ai = 1:NI
% % % %                 b = [b; ask_price_plus_total{ti,1}{ai,1}];
% % % %             end
% % % %             a = ask_price_plus_total{ti,1}{i,1}; ask_price_plus_total{ti,1}{i,1} = (a-refer_use);
%             
% %             b = [];
% %             for ai = 1:NI
% %                 b = [b; bid_price_plus_total{ti,1}{ai,1}];
% %             end
% %             b = b / sqrt(NI);
% %             a = bid_price_plus_total{ti,1}{i,1}; bid_price_plus_total{ti,1}{i,1} = (a-a(t+1)) + 100;
% %      
% %             b = [];
% %             for ai = 1:NI
% %                 b = [b; ask_price_plus_total{ti,1}{ai,1}];
% %             end
% %             b = b / sqrt(NI);
% %             a = ask_price_plus_total{ti,1}{i,1}; ask_price_plus_total{ti,1}{i,1} = (a-c(t+1)) + 100;
%             
%             
%             b = [];
%             for ai = 1:NI
%                 b = [b; spread_minus_total{ti,1}{ai,1}];
%             end
% %             b = b / sqrt(NI);
% %             a = spread_minus_total{ti,1}{i,1}; spread_minus_total{ti,1}{i,1} = (a-a(t+1)) + 100;
%             a = spread_minus_total{ti,1}{i,1}; spread_minus_total{ti,1}{i,1} = a;
%             
%             
%             b = [];
%             for ai = 1:NI
%                 b = [b; spread_plus_total{ti,1}{ai,1}];
%             end
% %             b = b / sqrt(NI);
% %             a = spread_plus_total{ti,1}{i,1}; spread_plus_total{ti,1}{i,1} = (a-a(t+1)) + 100;
%             a = spread_plus_total{ti,1}{i,1}; spread_plus_total{ti,1}{i,1} = a;
%             
% 
%             clear a
%         end
%     end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % delete non-value    
%     for ti=1:length(tick_bin_valid_save)
%         for i=1:NI
%             if sum(isnan(bid_price_minus_total{ti,1}{i,1})==1 | isnan(bid_price_plus_total{ti,1}{i,1})==1)>0
%                 for ii=1:5
%                     bid_vol_minus_total{ti,1}{i,ii} = 0;
%                     bid_vol_plus_total{ti,1}{i,ii} = 0;
%                     bid_gap_minus_total{ti,1}{i,ii} = 0;
%                     bid_gap_plus_total{ti,1}{i,ii} = 0;
%                     ask_vol_minus_total{ti,1}{i,ii} = 0;
%                     ask_vol_plus_total{ti,1}{i,ii} = 0;
%                     ask_gap_minus_total{ti,1}{i,ii} = 0;
%                     ask_gap_plus_total{ti,1}{i,ii} = 0;
%                 end
%                 for ii=1:4
%                     bid_vol_price_minus_avg_total{ti,1}{i,ii} = 0;
%                     bid_vol_price_plus_avg_total{ti,1}{i,ii} = 0;
%                     ask_vol_price_minus_avg_total{ti,1}{i,ii} = 0;
%                     ask_vol_price_plus_avg_total{ti,1}{i,ii} = 0;
%                 end
%                 for ii=1:18
%                     flow_minus_total{ti,1}{i,ii} = 0;
%                     flow_plus_total{ti,1}{i,ii} = 0;
%                 end
%                 
%                 bid_price_minus_total{ti,1}{i,1} = 0;
%                 bid_price_plus_total{ti,1}{i,1} = 0;
%                 ask_price_minus_total{ti,1}{i,1} = 0;
%                 ask_price_plus_total{ti,1}{i,1} = 0;
%                 spread_minus_total{ti,1}{i,1} = 0;
%                 spread_plus_total{ti,1}{i,1} = 0;
%                 update_time_gap_minus_total{ti,1}{i,1} = 0;
%                 update_time_gap_plus_total{ti,1}{i,1} = 0;
%             end
%         end
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    tick_bin_valid_save_go = zeros(length(tick_bin_valid_save),1);
    for ti=1:length(tick_bin_valid_save)
        for i=1:NI
            if sum(isnan(bid_price_minus_total{ti,1}{i,1})==1 | isnan(bid_price_plus_total{ti,1}{i,1})==1) > 0
                tick_bin_valid_save_go(ti,1) = tick_bin_valid_save_go(ti,1) + 1;
            end
        end
    end
    
    
    
    % needed to repair, just divide the values by # of tick_bin_valid
    % total aggregate
    
%     sum_minus = zeros(NI,1);    
%     sum_plus = zeros(NI,1);
    for ti=1:length(tick_bin_valid_save)
        for i=1:NI
            if bid_price_minus_num{ti,2}(i,1) > 0
                sum_minus(i,1) = sum_minus(i,1) + 1;
            end
            if bid_price_plus_num{ti,2}(i,1) > 0
                sum_plus(i,1) = sum_plus(i,1) + 1;
            end
        end
    end
    
    for ti=1:length(tick_bin_valid_save)
        if tick_bin_valid_save_go(ti) == 0
            for i=1:NI
                if sum_minus(i,1) > 0
                    for ii=1:5
                        bid_vol_minus_ta{1,1}{i,ii} = bid_vol_minus_ta{1,1}{i,ii} + bid_vol_minus_total{ti,1}{i,ii}; % / length(find(tick_bin_valid_save_go==0));
                        bid_gap_minus_ta{1,1}{i,ii} = bid_gap_minus_ta{1,1}{i,ii} + bid_gap_minus_total{ti,1}{i,ii};
                        ask_vol_minus_ta{1,1}{i,ii} = ask_vol_minus_ta{1,1}{i,ii} + ask_vol_minus_total{ti,1}{i,ii};
                        ask_gap_minus_ta{1,1}{i,ii} = ask_gap_minus_ta{1,1}{i,ii} + ask_gap_minus_total{ti,1}{i,ii};
                    end
%                     for ii=1:4
%                         bid_vol_price_minus_avg_ta{1,1}{i,ii} = bid_vol_price_minus_avg_ta{1,1}{i,ii} + bid_vol_price_minus_avg_total{ti,1}{i,ii};
%                         ask_vol_price_minus_avg_ta{1,1}{i,ii} = ask_vol_price_minus_avg_ta{1,1}{i,ii} + ask_vol_price_minus_avg_total{ti,1}{i,ii};
%                     end
%                     for ii=1:18
%                         flow_minus_ta{1,1}{i,ii} = flow_minus_ta{1,1}{i,ii} + flow_minus_total{ti,1}{i,ii};
%                     end

                    bid_price_minus_ta{1,1}{i,1} = bid_price_minus_ta{1,1}{i,1} + bid_price_minus_total{ti,1}{i,1};
                    ask_price_minus_ta{1,1}{i,1} = ask_price_minus_ta{1,1}{i,1} + ask_price_minus_total{ti,1}{i,1};
                    spread_minus_ta{1,1}{i,1} = spread_minus_ta{1,1}{i,1} + spread_minus_total{ti,1}{i,1};
%                     update_time_gap_minus_ta{1,1}{i,1} = update_time_gap_minus_ta{1,1}{i,1} + update_time_gap_minus_total{ti,1}{i,1};
                end
                if sum_plus(i,1) > 0
                    for ii=1:5
                        bid_vol_plus_ta{1,1}{i,ii} = bid_vol_plus_ta{1,1}{i,ii} + bid_vol_plus_total{ti,1}{i,ii};
                        bid_gap_plus_ta{1,1}{i,ii} = bid_gap_plus_ta{1,1}{i,ii} + bid_gap_plus_total{ti,1}{i,ii};
                        ask_vol_plus_ta{1,1}{i,ii} = ask_vol_plus_ta{1,1}{i,ii} + ask_vol_plus_total{ti,1}{i,ii};
                        ask_gap_plus_ta{1,1}{i,ii} = ask_gap_plus_ta{1,1}{i,ii} + ask_gap_plus_total{ti,1}{i,ii};
                    end
%                     for ii=1:4
%                         bid_vol_price_plus_avg_ta{1,1}{i,ii} = bid_vol_price_plus_avg_ta{1,1}{i,ii} + bid_vol_price_plus_avg_total{ti,1}{i,ii};
%                         ask_vol_price_plus_avg_ta{1,1}{i,ii} = ask_vol_price_plus_avg_ta{1,1}{i,ii} + ask_vol_price_plus_avg_total{ti,1}{i,ii};
%                     end
%                     for ii=1:18
%                         flow_plus_ta{1,1}{i,ii} = flow_plus_ta{1,1}{i,ii} + flow_plus_total{ti,1}{i,ii};
%                     end

                    bid_price_plus_ta{1,1}{i,1} = bid_price_plus_ta{1,1}{i,1} + bid_price_plus_total{ti,1}{i,1};
                    ask_price_plus_ta{1,1}{i,1} = ask_price_plus_ta{1,1}{i,1} + ask_price_plus_total{ti,1}{i,1};
                    spread_plus_ta{1,1}{i,1} = spread_plus_ta{1,1}{i,1} + spread_plus_total{ti,1}{i,1};
%                     update_time_gap_plus_ta{1,1}{i,1} = update_time_gap_plus_ta{1,1}{i,1} + update_time_gap_plus_total{ti,1}{i,1};
                end
            end
        end
    end
    
    for i=1:NI
        if sum_minus(i,1) > 0
            sum_minus(i,1) = 1;
        end
        if sum_plus(i,1) > 0
            sum_plus(i,1) = 1;
        end
    end
    
    sum_minus_ta(:,firm_ind) = sum_minus;
    sum_plus_ta(:,firm_ind) = sum_plus;
    
    firm_ind = firm_ind + 1;
    
    clear bid_vol_price_minus_avg_total bid_vol_price_plus_avg_total ask_vol_price_minus_avg_total ask_vol_price_plus_avg_total
    clear bid_vol_price_minus_avg_num bid_vol_price_plus_avg_num ask_vol_price_minus_avg_num ask_vol_price_plus_avg_num
    clear bid_vol_minus_total bid_vol_minus_num bid_vol_plus_total bid_vol_plus_num bid_gap_minus_total bid_gap_minus_num 
    clear bid_gap_plus_total bid_gap_plus_num ask_vol_minus_total ask_vol_minus_num ask_vol_plus_total ask_vol_plus_num
    clear ask_gap_minus_total ask_gap_minus_num ask_gap_plus_total ask_gap_plus_num flow_minus_total flow_minus_num
    clear flow_plus_total flow_plus_num bid_price_minus_total bid_price_minus_num bid_price_plus_total bid_price_plus_num
    clear ask_price_minus_total ask_price_minus_num ask_price_plus_total ask_price_plus_num spread_minus_total spread_minus_num
    clear spread_plus_total spread_plus_num update_time_gap_minus_total update_time_gap_minus_num
    clear update_time_gap_plus_total update_time_gap_plus_num
    clear sum_minus sum_plus
end

for i=1:NI
    for ii=1:5
        bid_vol_minus_ta{1,1}{i,ii} = bid_vol_minus_ta{1,1}{i,ii} / num_event_minus(i,1);
        bid_vol_plus_ta{1,1}{i,ii} = bid_vol_plus_ta{1,1}{i,ii} / num_event_plus(i,1);
%         bid_vol_plus_ta{1,1}{i,ii} = bid_vol_plus_ta{1,1}{i,ii} / (firm_end-firm_start+1);
        bid_gap_minus_ta{1,1}{i,ii} = bid_gap_minus_ta{1,1}{i,ii} / num_event_minus(i,1);
        bid_gap_plus_ta{1,1}{i,ii} = bid_gap_plus_ta{1,1}{i,ii} / num_event_plus(i,1);
        ask_vol_minus_ta{1,1}{i,ii} = ask_vol_minus_ta{1,1}{i,ii}/ num_event_minus(i,1);
        ask_vol_plus_ta{1,1}{i,ii} = ask_vol_plus_ta{1,1}{i,ii} / num_event_plus(i,1);
        ask_gap_minus_ta{1,1}{i,ii} = ask_gap_minus_ta{1,1}{i,ii} / num_event_minus(i,1);
        ask_gap_plus_ta{1,1}{i,ii} = ask_gap_plus_ta{1,1}{i,ii} / num_event_plus(i,1);


%         bid_vol_minus_ta{1,1}{i,ii} = bid_vol_minus_ta{1,1}{i,ii} / (firm_end-firm_start+1);
%         bid_vol_plus_ta{1,1}{i,ii} = bid_vol_plus_ta{1,1}{i,ii} / (firm_end-firm_start+1);
% %         bid_vol_plus_ta{1,1}{i,ii} = bid_vol_plus_ta{1,1}{i,ii} / (firm_end-firm_start+1);
%         bid_gap_minus_ta{1,1}{i,ii} = bid_gap_minus_ta{1,1}{i,ii} / (firm_end-firm_start+1);
%         bid_gap_plus_ta{1,1}{i,ii} = bid_gap_plus_ta{1,1}{i,ii} / (firm_end-firm_start+1);
%         ask_vol_minus_ta{1,1}{i,ii} = ask_vol_minus_ta{1,1}{i,ii} / (firm_end-firm_start+1);
%         ask_vol_plus_ta{1,1}{i,ii} = ask_vol_plus_ta{1,1}{i,ii} / (firm_end-firm_start+1);
%         ask_gap_minus_ta{1,1}{i,ii} = ask_gap_minus_ta{1,1}{i,ii} / (firm_end-firm_start+1);
%         ask_gap_plus_ta{1,1}{i,ii} = ask_gap_plus_ta{1,1}{i,ii} / (firm_end-firm_start+1);
    end
%     for ii=1:4
%         bid_vol_price_minus_avg_ta{1,1}{i,ii} = bid_vol_price_minus_avg_ta{1,1}{i,ii} / num_event_minus(i,1);
%         bid_vol_price_plus_avg_ta{1,1}{i,ii} = bid_vol_price_plus_avg_ta{1,1}{i,ii} / num_event_plus(i,1);
%         ask_vol_price_minus_avg_ta{1,1}{i,ii} = ask_vol_price_minus_avg_ta{1,1}{i,ii} / num_event_minus(i,1);
%         ask_vol_price_plus_avg_ta{1,1}{i,ii} = ask_vol_price_plus_avg_ta{1,1}{i,ii} / num_event_plus(i,1);
%     end
%     for ii=1:18
%         flow_minus_ta{1,1}{i,ii} = flow_minus_ta{1,1}{i,ii} / num_event_minus(i,1);
%         flow_plus_ta{1,1}{i,ii} = flow_plus_ta{1,1}{i,ii}/ num_event_plus(i,1);
%     end
    bid_price_minus_ta{1,1}{i,1} = bid_price_minus_ta{1,1}{i,1} / num_event_minus(i,1);
    bid_price_plus_ta{1,1}{i,1} = bid_price_plus_ta{1,1}{i,1} / num_event_plus(i,1);
    ask_price_minus_ta{1,1}{i,1} = ask_price_minus_ta{1,1}{i,1} / num_event_minus(i,1);
    ask_price_plus_ta{1,1}{i,1} = ask_price_plus_ta{1,1}{i,1} / num_event_plus(i,1);
%     bid_price_minus_ta{1,1}{i,1} = bid_price_minus_ta{1,1}{i,1} / sum(sum_minus_ta(i,:)) / num_event_minus(i,1);
%     bid_price_plus_ta{1,1}{i,1} = bid_price_plus_ta{1,1}{i,1} / sum(sum_plus_ta(i,:)) / num_event_plus(i,1);
%     ask_price_minus_ta{1,1}{i,1} = ask_price_minus_ta{1,1}{i,1} / sum(sum_minus_ta(i,:)) / num_event_minus(i,1);
%     ask_price_plus_ta{1,1}{i,1} = ask_price_plus_ta{1,1}{i,1} / sum(sum_plus_ta(i,:)) / num_event_plus(i,1);
    spread_minus_ta{1,1}{i,1} = spread_minus_ta{1,1}{i,1} / num_event_minus(i,1);
    spread_plus_ta{1,1}{i,1} = spread_plus_ta{1,1}{i,1} / num_event_plus(i,1);
%     update_time_gap_minus_ta{1,1}{i,1} = update_time_gap_minus_ta{1,1}{i,1} / num_event_minus(i,1);
%     update_time_gap_plus_ta{1,1}{i,1} = update_time_gap_plus_ta{1,1}{i,1} / num_event_plus(i,1);
end


% save data
save(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/result_AB_Burst/z15o_basket_%d_%d_%d_%d.mat',firm_start,firm_end,day_start,day_end))

% % % %%
% % % load('/home/minyoung/data_minyoung/Research/Price_Impact/result_AB_Burst/z15o.mat')

%%
% % % %% load data
% % % clear;clc;
% % % load('/home/minyoung/data_minyoung/Research/Price_Impact/ver6_repair_g_to_script_ini.mat')
% % % % load('/home/minyoung/data_minyoung/Research/Price_Impact/ver6_repair_g_to_script.mat')
% % % % load('/home/minyoung/data_minyoung/Research/Price_Impact/ver6_repair_g_to_1000.mat')
% % % 
% % % %% # of event
% % % k=cell(NI,1);
% % % for i=1:NI
% % %     k{i,1}=zeros(1,510);
% % % end
% % % tick_bin_valid_save=[];
% % % for firm=firm_start:firm_end
% % %     for day=day_start:day_end
% % %         if exist(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_repair_g_k_%d_%d_%d.mat',com_num(firm),day,id)) && day_ind ~=26
% % %             load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_repair_g_k_%d_%d_%d.mat',com_num(firm),day,id))
% % %             for i=1:NI
% % %                 k{i,1} = k{i,1} + k1{i,1} + k2{i,1};
% % %             end
% % %         end
% % %         day_ind = day_ind + 1;
% % %     end
% % % end
% % % 
% % % % k=k(rearrange_ind);
% % % 
% % % xtime = linspace(8,16.5,510);
% % % figure
% % % set(gcf,'color','w')
% % % for i=1:NI
% % %     subplot(3,5,i)
% % %     plot(xtime,k{i,1},'k')
% % %     xlim([8 16.5])
% % %     s(i)=sum(k{i,1});
% % % end
% % % figure
% % % set(gcf,'color','w')
% % % plot(s/(len_day),'-*','LineWidth',2)
% % % hold on
% % % plot([1,2,4,7,11],s([1,2,4,7,11])/(len_day),'ko--','MarkerSize',8,'MarkerFaceColor',[0 0 0])
% % % % scatter([1,2,4,7,11],s([1,2,4,7,11])/(day_end-day_start+1),'ko','filled')
% % % set(gca,'yscale','log')
% % % xlabel('Type of event')
% % % ylabel('Number of event')
% % % legend('All event','Type A')
% % % 
% % % 
% % % %% spread (1 figure, normalization x)
% % % iter=iter_set;
% % % figure100 = figure('Position',[1 1 1000 800],'Color',[1 1 1]);
% % % set(gcf,'color','w')
% % % % set(f1,'Visible','off')
% % % % set(f1, 'PaperUnits', 'inches');
% % % % x_width=14 ;y_width=7;
% % % % set(f1, 'PaperPosition', [0 0 x_width y_width]);
% % % for i=1:NI
% % %     if i==1
% % %         ind = 1;
% % %     elseif i>1 && i<4
% % %         ind = 2;
% % %         col = (i-1)/3;
% % %     elseif i>3 && i<7
% % %         ind = 3;
% % %         col = (i-3)/4;
% % %     elseif i>6 && i<11
% % %         ind = 4;
% % %         col = (i-6)/5;
% % %     elseif i>10
% % %         ind = 5;
% % %         col = (i-10)/6;
% % %     end
% % %     
% % %     if i==1 || i==2 || i==4 || i==7 || i==11
% % %         flag_g = 1;
% % %     else
% % %         flag_g = 2;
% % %     end
% % %     subplot(2,3,ind)
% % % %         hold on
% % % %         plot(0,spread_plus_total{iter,1}{i,1}(t+1)./spread_plus_num{iter,1}(i,1),'k*','MarkerSize',10)
% % % %         hold on
% % %     if flag_g==1
% % % %             plot([-t:1:-1,1:t],spread_plus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_plus_num{iter,1}(i,1),'k.')
% % %         plot(-t:t,spread_plus_ta{iter,1}{i,1}(1:2*t+1),'k-') %  - spread_plus_ta{iter,1}{i,1}(t+2) + 8
% % %         hold on
% % %         plot(0,spread_plus_ta{iter,1}{i,1}(t+1),'k*','MarkerSize',10)
% % %         hold on
% % %         plot(0,spread_plus_ta{iter,1}{i,1}(t+2),'k*','MarkerSize',10)
% % %     elseif flag_g==2
% % % %             plot([-t:1:-1,1:t],spread_plus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_plus_num{iter,1}(i,1),'b')
% % %         plot(-t:t,spread_plus_ta{iter,1}{i,1}(1:2*t+1),'Color',[1 col col]) %  - spread_plus_ta{iter,1}{i,1}(t+2) + 8
% % %         hold on
% % %         plot(0,spread_plus_ta{iter,1}{i,1}(t+1),'o','Color',[1 col col],'MarkerSize',10)
% % %         hold on
% % %         plot(0,spread_plus_ta{iter,1}{i,1}(t+2),'o','Color',[1 col col],'MarkerSize',10)
% % %     end
% % %         ylim([1.3 6.4])
% % % %         ylim([1 9])
% % % %         set(gca,'xscale','log')
% % % %         set(gca,'yscale','log')
% % % %         hold on
% % % %     hold on
% % % %     plot([-t:1:-1,1:t],spread_minus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_minus_num{iter,1}(i,1),'b.')
% % % %     legend('minus initiated market order','plus initiated market order')
% % %     xlabel('\tau')
% % %     if ind==1
% % %         ylabel('G_s(\tau|\Delta=1)')
% % %     elseif ind==2
% % %         ylabel('G_s(\tau|\Delta=2)')
% % %     elseif ind==3
% % %         ylabel('G_s(\tau|\Delta=3)')
% % %     elseif ind==4
% % %         ylabel('G_s(\tau|\Delta=4)')
% % %     elseif ind==5
% % %         ylabel('G_s(\tau|\Delta=5)')
% % %     end
% % %     
% % %     if i==1
% % %         text(-0.1,1.1,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
% % %     elseif i==2
% % %         text(-0.1,1.1,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
% % %     elseif i==4
% % %         text(-0.1,1.1,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top')
% % %     elseif i==7
% % %         text(-0.1,1.1,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top')
% % %     elseif i==11
% % %         text(-0.1,1.1,'(e)','Units', 'Normalized', 'VerticalAlignment', 'Top')
% % %     end
% % % end
% % % 
% % % dx=0.05;
% % % dy=0.05;
% % % x=(1-5*dx)/2.9;
% % % y=(1-5*dy)/2;
% % % % dxx=0.01;
% % % % dyy=0.01;
% % % AxesHandle=findobj(figure100,'Type','axes');
% % % set(AxesHandle(5),'Position',[0.01+dx,y+4*dy,x,y]);
% % % set(AxesHandle(4),'Position',[0.01+x+2.5*dx,y+4*dy,x,y]);
% % % set(AxesHandle(3),'Position',[0.01+2*x+4*dx,y+4*dy,x,y]);
% % % set(AxesHandle(2),'Position',[0.01+dx,dy,x,y]);
% % % set(AxesHandle(1),'Position',[0.01+x+2.5*dx,dy,x,y]);
% % % 
% % % % iter=iter_set;
% % % % f2=figure;
% % % % set(gcf,'color','w')
% % % % % set(f1,'Visible','off')
% % % % % set(f1, 'PaperUnits', 'inches');
% % % % % x_width=14 ;y_width=7;
% % % % % set(f1, 'PaperPosition', [0 0 x_width y_width]);
% % % % for i=1:NI
% % % %     if i>5 && i<10
% % % %         ind = i-4;
% % % %         flag = 1;
% % % %         flag_g = 1;
% % % %     elseif i>9
% % % %         ind = i-8;
% % % %         flag = 1;
% % % %         flag_g = 2;
% % % %     elseif i==1
% % % %         ind = 1;
% % % %         flag = 1;
% % % %         flag_g = 1;
% % % %     elseif i>1 && i<6
% % % %         flag = 0;
% % % %     end
% % % %     if flag==1
% % % %         subplot(1,5,ind)
% % % % %         plot(0,spread_plus_total{iter,1}{i,1}(t+1)./spread_plus_num{iter,1}(i,1),'k*','MarkerSize',10)
% % % % %         hold on
% % % %         if flag_g==1
% % % % %             plot([-t:1:-1,1:t],spread_plus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_plus_num{iter,1}(i,1),'k.')
% % % %             plot(-t:t,spread_plus_ta{iter,1}{i,1}(1:2*t+1),'k-') %  - spread_plus_ta{iter,1}{i,1}(t+2) + 8
% % % %             hold on
% % % %             plot(0,spread_plus_ta{iter,1}{i,1}(t+1),'k*','MarkerSize',10)
% % % %             hold on
% % % %             plot(0,spread_plus_ta{iter,1}{i,1}(t+2),'k*','MarkerSize',10)
% % % %         elseif flag_g==2
% % % % %             plot([-t:1:-1,1:t],spread_plus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_plus_num{iter,1}(i,1),'b')
% % % %             plot(-t:t,spread_plus_ta{iter,1}{i,1}(1:2*t+1),'b-') %  - spread_plus_ta{iter,1}{i,1}(t+2) + 8
% % % %             hold on
% % % %             plot(0,spread_plus_ta{iter,1}{i,1}(t+1),'bo','MarkerSize',10)
% % % %             hold on
% % % %             plot(0,spread_plus_ta{iter,1}{i,1}(t+2),'bo','MarkerSize',10)
% % % %         end
% % % % %         ylim([1.3 8])
% % % % %         ylim([99.5 106])
% % % % %         set(gca,'xscale','log')
% % % % %         set(gca,'yscale','log')
% % % % %         grid on
% % % % %         hold on
% % % %     %     hold on
% % % %     %     plot([-t:1:-1,1:t],spread_minus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_minus_num{iter,1}(i,1),'b.')
% % % %     %     legend('minus initiated market order','plus initiated market order')
% % % %         xlabel('lag')
% % % %         ylabel('Spread')
% % % %     end
% % % % end
% % % 
% % % 
% % % 
% % % 
% % % 
% % % %% spread (1 figure - considering normal spread)
% % % iter=iter_set;
% % % f1=figure;
% % % set(gcf,'color','w')
% % % % set(f1,'Visible','off')
% % % % set(f1, 'PaperUnits', 'inches');
% % % % x_width=14 ;y_width=7;
% % % % set(f1, 'PaperPosition', [0 0 x_width y_width]);
% % % for i=1:NI
% % %     if i==1
% % %         ind = 1;
% % %     elseif i>1 && i<4
% % %         ind = 2;
% % %         col = (i-1)/3;
% % %     elseif i>3 && i<7
% % %         ind = 3;
% % %         col = (i-3)/4;
% % %     elseif i>6 && i<11
% % %         ind = 4;
% % %         col = (i-6)/5;
% % %     elseif i>10
% % %         ind = 5;
% % %         col = (i-10)/6;
% % %     end
% % %     
% % %     if i==1 || i==2 || i==4 || i==7 || i==11
% % %         flag_g = 1;
% % %     else
% % %         flag_g = 2;
% % %     end
% % %     
% % %     subplot(1,5,ind)
% % % %         plot(0,spread_plus_total{iter,1}{i,1}(t+1)./spread_plus_num{iter,1}(i,1),'k*','MarkerSize',10)
% % % %         hold on
% % %     if flag_g==1
% % % %             plot([-t:1:-1,1:t],spread_plus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_plus_num{iter,1}(i,1),'k.')
% % %         plot(-t:t,spread_plus_ta{iter,1}{i,1}(1:2*t+1) - mean(spread_plus_ta{iter,1}{i,1}(1:20)),'k-') %  - spread_plus_ta{iter,1}{i,1}(t+2) + 8
% % %         hold on
% % %         plot(0,spread_plus_ta{iter,1}{i,1}(t+1) - mean(spread_plus_ta{iter,1}{i,1}(1:20)),'k*','MarkerSize',10)
% % %         hold on
% % %         plot(0,spread_plus_ta{iter,1}{i,1}(t+2) - mean(spread_plus_ta{iter,1}{i,1}(1:20)),'k*','MarkerSize',10)
% % %     elseif flag_g==2
% % % %             plot([-t:1:-1,1:t],spread_plus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_plus_num{iter,1}(i,1),'b')
% % %         plot(-t:t,spread_plus_ta{iter,1}{i,1}(1:2*t+1) - mean(spread_plus_ta{iter,1}{i,1}(1:20)),'Color',[1 col col]) %  - spread_plus_ta{iter,1}{i,1}(t+2) + 8
% % %         hold on
% % %         plot(0,spread_plus_ta{iter,1}{i,1}(t+1) - mean(spread_plus_ta{iter,1}{i,1}(1:20)),'o','Color',[1 col col],'MarkerSize',10)
% % %         hold on
% % %         plot(0,spread_plus_ta{iter,1}{i,1}(t+2) - mean(spread_plus_ta{iter,1}{i,1}(1:20)),'o','Color',[1 col col],'MarkerSize',10)
% % %     end
% % % %         ylim([99.5 106])
% % % %         ylim([1 9])
% % %         set(gca,'xscale','log')
% % %         set(gca,'yscale','log')
% % % %         hold on
% % % %     hold on
% % % %     plot([-t:1:-1,1:t],spread_minus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_minus_num{iter,1}(i,1),'b.')
% % % %     legend('minus initiated market order','plus initiated market order')
% % %     xlabel('lag')
% % %     ylabel('Spread')
% % % end
% % % 
% % % 
% % % % iter=iter_set;
% % % % f2=figure;
% % % % set(gcf,'color','w')
% % % % % set(f1,'Visible','off')
% % % % % set(f1, 'PaperUnits', 'inches');
% % % % % x_width=14 ;y_width=7;
% % % % % set(f1, 'PaperPosition', [0 0 x_width y_width]);
% % % % for i=1:NI
% % % %     if i>5 && i<10
% % % %         ind = i-4;
% % % %         flag = 1;
% % % %         flag_g = 1;
% % % %     elseif i>9
% % % %         ind = i-8;
% % % %         flag = 1;
% % % %         flag_g = 2;
% % % %     elseif i==1
% % % %         ind = 1;
% % % %         flag = 1;
% % % %         flag_g = 1;
% % % %     elseif i>1 && i<6
% % % %         flag = 0;
% % % %     end
% % % %     if flag==1
% % % %         subplot(1,5,ind)
% % % % %         plot(0,spread_plus_total{iter,1}{i,1}(t+1)./spread_plus_num{iter,1}(i,1),'k*','MarkerSize',10)
% % % % %         hold on
% % % %         if flag_g==1
% % % % %             plot([-t:1:-1,1:t],spread_plus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_plus_num{iter,1}(i,1),'k.')
% % % %             plot(-t:t,spread_minus_ta{iter,1}{i,1}(1:2*t+1) - mean(spread_minus_ta{iter,1}{i,1}(1:20)),'k-') %  - spread_plus_ta{iter,1}{i,1}(t+2) + 8
% % % %             hold on
% % % %             plot(0,spread_minus_ta{iter,1}{i,1}(t+1) - mean(spread_minus_ta{iter,1}{i,1}(1:20)),'k*','MarkerSize',10)
% % % %             hold on
% % % %             plot(0,spread_minus_ta{iter,1}{i,1}(t+2) - mean(spread_minus_ta{iter,1}{i,1}(1:20)),'k*','MarkerSize',10)
% % % %         elseif flag_g==2
% % % % %             plot([-t:1:-1,1:t],spread_plus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_plus_num{iter,1}(i,1),'b')
% % % %             plot(-t:t,spread_minus_ta{iter,1}{i,1}(1:2*t+1) - mean(spread_minus_ta{iter,1}{i,1}(1:20)),'b-') %  - spread_plus_ta{iter,1}{i,1}(t+2) + 8
% % % %             hold on
% % % %             plot(0,spread_minus_ta{iter,1}{i,1}(t+1) - mean(spread_minus_ta{iter,1}{i,1}(1:20)),'bo','MarkerSize',10)
% % % %             hold on
% % % %             plot(0,spread_minus_ta{iter,1}{i,1}(t+2) - mean(spread_minus_ta{iter,1}{i,1}(1:20)),'bo','MarkerSize',10)
% % % %         end
% % % % %         ylim([1.3 8])
% % % % %         ylim([99.5 106])
% % % %         set(gca,'xscale','log')
% % % %         set(gca,'yscale','log')
% % % % %         grid on
% % % % %         hold on
% % % %     %     hold on
% % % %     %     plot([-t:1:-1,1:t],spread_minus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_minus_num{iter,1}(i,1),'b.')
% % % %     %     legend('minus initiated market order','plus initiated market order')
% % % %         xlabel('lag')
% % % %         ylabel('Spread')
% % % %     end
% % % % end
% % % 
% % % %% time interval (1 figure)
% % % iter=iter_set;
% % % f1=figure;
% % % set(gcf,'color','w')
% % % % set(f1,'Visible','off')
% % % % set(f1, 'PaperUnits', 'inches');
% % % % x_width=14 ;y_width=7;
% % % % set(f1, 'PaperPosition', [0 0 x_width y_width]);
% % % for i=1:NI
% % %     if i==1
% % %         ind = 1;
% % %     elseif i>1 && i<4
% % %         ind = 2;
% % %         col = (i-1)/3;
% % %     elseif i>3 && i<7
% % %         ind = 3;
% % %         col = (i-3)/4;
% % %     elseif i>6 && i<11
% % %         ind = 4;
% % %         col = (i-6)/5;
% % %     elseif i>10
% % %         ind = 5;
% % %         col = (i-10)/6;
% % %     end
% % %     
% % %     if i==1 || i==2 || i==4 || i==7 || i==11
% % %         flag_g = 1;
% % %     else
% % %         flag_g = 2;
% % %     end
% % %     subplot(1,5,ind)
% % % %         plot(0,spread_plus_total{iter,1}{i,1}(t+1)./spread_plus_num{iter,1}(i,1),'k*','MarkerSize',10)
% % % %         hold on
% % %     if flag_g==1
% % % %             plot([-t:1:-1,1:t],spread_plus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_plus_num{iter,1}(i,1),'k.')
% % %         plot(-t:t,update_time_gap_plus_ta{iter,1}{i,1}(1:2*t+1),'k-') %  - spread_plus_ta{iter,1}{i,1}(t+2) + 8
% % %         hold on
% % %         plot(0,update_time_gap_plus_ta{iter,1}{i,1}(t+1),'k*','MarkerSize',10)
% % %         hold on
% % %         plot(0,update_time_gap_plus_ta{iter,1}{i,1}(t+2),'k*','MarkerSize',10)
% % %     elseif flag_g==2
% % % %             plot([-t:1:-1,1:t],spread_plus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_plus_num{iter,1}(i,1),'b')
% % %         plot(-t:t,update_time_gap_plus_ta{iter,1}{i,1}(1:2*t+1),'Color',[1 col col]) %  - spread_plus_ta{iter,1}{i,1}(t+2) + 8
% % %         hold on
% % %         plot(0,update_time_gap_plus_ta{iter,1}{i,1}(t+1),'o','Color',[1 col col],'MarkerSize',10)
% % %         hold on
% % %         plot(0,update_time_gap_plus_ta{iter,1}{i,1}(t+2),'o','Color',[1 col col],'MarkerSize',10)
% % %     end
% % %         xlim([-100 100])
% % %         ylim([0 1700])
% % % %         set(gca,'xscale','log')
% % % %         set(gca,'yscale','log')
% % % %         hold on
% % % %     hold on
% % % %     plot([-t:1:-1,1:t],spread_minus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_minus_num{iter,1}(i,1),'b.')
% % % %     legend('minus initiated market order','plus initiated market order')
% % %     xlabel('lag')
% % %     ylabel('Spread')
% % % end
% % % 
% % % 

%%
% clear;clc;

firm_start = 1;
firm_end = 65;%70;
day_start = 32;
day_end = 61;

addpath('/home/minyoung/data_minyoung/Research/Matlab_function/linspecer')
load(sprintf('./result_AB_Burst/z15o_basket_%d_%d_%d_%d.mat',firm_start,firm_end,day_start,day_end));
% load(sprintf('./result_AB_Burst/z15o_%d_%d_%d_%d.mat',firm_start,firm_end,day_start,day_end));

ci = 2;
lind = 1;
lw = 1;
delta_t = 100;
we = t+1+delta_t;

C2 = linspecer(3);
C4 = linspecer(3);
C2(1,:) = C2(3,:);
C3 = linspecer(10);

c1 = 0;
c2 = 0;
c3 = 0;

% c1 = C2(3,1);
% c2 = C2(3,2);
% c3 = C2(3,3);

% C2(1,:) = [1, 0, 0];
% C2(2,:) = [0, 0, 1];
% C2 = linspecer(10);
% C2(1,:) = C2(3,:);


%% 

addpath('/home/minyoung/data_minyoung/Research/Matlab_function/linspecer')
% C = linspecer(4);
C(1,:) = [0,0,1];
C(2,:) = [1,0,0];
C = [57,106,177;204,37,41]/255;
iter=iter_set;
% figure100 = figure('Position',[-5 1 1000 800],'Color',[1 1 1]);
figure1 = figure(1);
set(gcf,'color','w')
% set(f1,'Visible','off')
% set(f1, 'PaperUnits', 'inches');
% x_width=14 ;y_width=7;
% set(f1, 'PaperPosition', [0 0 x_width y_width]);


for i=1:NI-2
    if i==1
        ind = 2;
    elseif i>1 && i<4
        ind = 3;
%         col = (i-1)/3;
    elseif i>3 && i<6
        ind = 4;
%         col = (i-3)/4;
    elseif i>5 && i<8
        ind = 5;
%         col = (i-6)/5;
    elseif i>7 && i<10
        ind = 6;
%         col = (i-10)/6;
    elseif i == 10
        ind = 1;
%         col = (i-10)/6;
    end
    
    if i==2 || i==4 || i==6 || i==8
        flag_g = 1;
    elseif i==3 || i==5 || i==7 || i==9
        flag_g = 2;
    elseif i==10
        flag_g = 0;
    elseif i == 1
        flag_g = 3;
    else
        flag_g = 4;
    end
    
    refer_use =  ask_price_plus_ta{iter,1}{i,1}(t+1);
%     refer_use = 0;
   

    subplot(2,3,ind)

    h10 = plot([-t,t],[0,0],':','color',[1/2 1/2 1/2]);
    hold on
    
%     co = [0.0336;0.0380;0.6581;0.0493;0.8230;0.0627;0.9409;0.0772;1.0346;0];
    co = [0.0305;0.0335;0.5364;0.0442;0.6474;0.0558;0.7315;0.0679;0.7851;0];
    
    if flag_g == 1
        t1 = plot(1,co(i,1),'s','color',[0 0 0],'LineWidth',1,'MarkerSize',10);
    elseif flag_g == 2
        t2 = plot(1,co(i,1),'x','color',[0 0 0],'LineWidth',1,'MarkerSize',10);
    elseif flag_g == 0
        t3 = plot(1,co(i,1),'o','color',[0 0 0],'LineWidth',1,'MarkerSize',10);
    elseif flag_g == 3
        t4 = plot(1,co(i,1),'d','color',[0 0 0],'LineWidth',1,'MarkerSize',10);
    end
    hold on

    if flag_g==1 % type A
        h1 = plot(-t:t,bid_price_plus_ta{iter,1}{i,1} - refer_use,'-','color',C(2,:),'LineWidth',1);
        bb = bid_price_plus_ta{iter,1}{i,1} - refer_use + 3;
        hold on
        h2 = plot(-t:t,ask_price_plus_ta{iter,1}{i,1} - refer_use,'-','color',C(1,:),'LineWidth',1);
        aa = ask_price_plus_ta{iter,1}{i,1} - refer_use + 3;
        hold on
%         plot(0,bid_price_plus_ta{iter,1}{i,1}(t+1) - refer_use,'*','color',C(2,:),'MarkerSize',10)
%         hold on
%         plot(0,ask_price_plus_ta{iter,1}{i,1}(t+1) - refer_use,'*','color',C(1,:),'MarkerSize',10)
%         hold on
%         plot(1,bid_price_plus_ta{iter,1}{i,1}(t+2) - refer_use,'*','color',C(2,:),'MarkerSize',10)
%         hold on
%         plot(1,ask_price_plus_ta{iter,1}{i,1}(t+2) - refer_use,'*','color',C(1,:),'MarkerSize',10)
    elseif flag_g==2 % type B
        h3 = plot(-t:t,bid_price_plus_ta{iter,1}{i,1} - refer_use,'--','color',C(2,:),'LineWidth',1);
        bb = bid_price_plus_ta{iter,1}{i,1} - refer_use + 3;
        hold on
        h4 = plot(-t:t,ask_price_plus_ta{iter,1}{i,1} - refer_use,'--','color',C(1,:),'LineWidth',1);
        aa = ask_price_plus_ta{iter,1}{i,1} - refer_use + 3;
        hold on
%         plot(-t:t,bid_price_plus_ta{iter,1}{i,1} - refer_use,'r.','LineWidth',1)
%         hold on
%         plot(-t:t,ask_price_plus_ta{iter,1}{i,1} - refer_use,'b.','LineWidth',1)
%         hold on
%         plot(0,bid_price_plus_ta{iter,1}{i,1}(t+1) - refer_use,'o','color',C(2,:),'MarkerSize',10)
%         hold on
%         plot(0,ask_price_plus_ta{iter,1}{i,1}(t+1) - refer_use,'o','color',C(1,:),'MarkerSize',10)
%         hold on
%         plot(1,bid_price_plus_ta{iter,1}{i,1}(t+2) - refer_use,'o','color',C(2,:),'MarkerSize',10)
%         hold on
%         plot(1,ask_price_plus_ta{iter,1}{i,1}(t+2) - refer_use,'o','color',C(1,:),'MarkerSize',10)
    elseif flag_g == 0% type Zero
        h5 = plot(-t:t,bid_price_plus_ta{iter,1}{i,1} - refer_use,'-.','color',C(2,:),'LineWidth',1);
        bb = bid_price_plus_ta{iter,1}{i,1} - refer_use + 3;
        hold on
        h6 = plot(-t:t,ask_price_plus_ta{iter,1}{i,1} - refer_use,'-.','color',C(1,:),'LineWidth',1);
        aa = ask_price_plus_ta{iter,1}{i,1} - refer_use + 3;
        hold on
    elseif flag_g == 3% type One
        h7 = plot(-t:t,bid_price_plus_ta{iter,1}{i,1} - refer_use,'-.','color',C(2,:),'LineWidth',1);
        bb = bid_price_plus_ta{iter,1}{i,1} - refer_use + 3;
        hold on
        h8 = plot(-t:t,ask_price_plus_ta{iter,1}{i,1} - refer_use,'-.','color',C(1,:),'LineWidth',1);
        aa = ask_price_plus_ta{iter,1}{i,1} - refer_use + 3;
        hold on
%         plot(0,bid_price_plus_ta{iter,1}{i,1}(t+1) - refer_use,'^','color',C(2,:),'MarkerSize',10)
%         hold on
%         plot(0,ask_price_plus_ta{iter,1}{i,1}(t+1) - refer_use,'^','color',C(1,:),'MarkerSize',10)
%         hold on
%         plot(1,bid_price_plus_ta{iter,1}{i,1}(t+2) - refer_use,'^','color',C(2,:),'MarkerSize',10)
%         hold on
%         plot(1,ask_price_plus_ta{iter,1}{i,1}(t+2) - refer_use,'^','color',C(1,:),'MarkerSize',10)
    end
    ylim([-3 5])
    xlim([-delta_t delta_t])
    ra_aa(i) = mean(aa(end-49:end)) - mean(aa(1:50));
    ra_bb(i) = mean(bb(end-49:end)) - mean(bb(1:50));
%     ylim([-2.5 5.2])
%     xlim([-.5 100])
%         set(gca,'xscale','log')
%         set(gca,'yscale','log')
%         hold on
%     hold on
%     plot([-t:1:-1,1:t],spread_minus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_minus_num{iter,1}(i,1),'b.')
%     legend('minus initiated market order','plus initiated market order')
    xlabel('\tau')
    if i == 1 || i == 3 || i== 5 || i == 7 || i == 9 || i == 10
        if ind==2
            ylabel('G_a(\tau|\Deltaa_0=1),  G_b(\tau|\Deltaa_0=1)')
            lh = legend([h7 h8],{'Bid (Type One)','Ask (Type One)'});
            legend boxoff
            set(lh, 'Position', [.197 .92 .5 .01])
        elseif ind==3
            ylabel('G_a(\tau|\Deltaa_0=2),  G_b(\tau|\Deltaa_0=2)')
            lh = legend([h1 h2 h3 h4],{'Bid (Type A)','Ask (Type A)','Bid (Type B)','Ask (Type B)'});
            legend boxoff
            set(lh, 'Position', [.53 .90 .5 .01])
        elseif ind==4
            ylabel('G_a(\tau|\Deltaa_0=3),  G_b(\tau|\Deltaa_0=3)')
            lh = legend([h1 h2 h3 h4],{'Bid (Type A)','Ask (Type A)','Bid (Type B)','Ask (Type B)'});
            legend boxoff
            set(lh, 'Position', [-.135 .41 .5 .01])
        elseif ind==5
            ylabel('G_a(\tau|\Deltaa_0=4),  G_b(\tau|\Deltaa_0=4)')
            lh = legend([h1 h2 h3 h4],{'Bid (Type A)','Ask (Type A)','Bid (Type B)','Ask (Type B)'});
            legend boxoff
            set(lh, 'Position', [.197 .41 .5 .01])
        elseif ind==6
            ylabel('G_a(\tau|\Deltaa_0=5),  G_b(\tau|\Deltaa_0=5)')
            lh = legend([h1 h2 h3 h4],{'Bid (Type A)','Ask (Type A)','Bid (Type B)','Ask (Type B)'});
            legend boxoff
            set(lh, 'Position', [.53 .41 .5 .01])
        elseif ind==1
            ylabel('G_a(\tau|\Deltaa_0=0),  G_b(\tau|\Deltaa_0=0)')
            lh = legend([h5 h6],{'Bid (Type Zero)','Ask (Type Zero)'});
            legend boxoff
            set(lh, 'Position', [-.128 .92 .5 .01])
        end
        if i==1
            text(-0.1,1.1,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
        elseif i==3
            text(-0.1,1.1,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top')
        elseif i==5
            text(-0.1,1.1,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top')
        elseif i==7
            text(-0.1,1.1,'(e)','Units', 'Normalized', 'VerticalAlignment', 'Top')
        elseif i==9
            text(-0.1,1.1,'(f)','Units', 'Normalized', 'VerticalAlignment', 'Top')
        elseif i==10
            text(-0.1,1.1,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
        end
    end
end
dx=0.05;
dy=0.05;
x=(1-5*dx)/2.9;
y=(1-5*dy)/1.9;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure1,'Type','axes');
set(AxesHandle(1),'Position',[dx,y+3.45*dy,x,y]);
set(AxesHandle(6),'Position',[x+2.5*dx,y+3.45*dy,x,y]);
set(AxesHandle(5),'Position',[2*x+4*dx,y+3.45*dy,x,y]);
set(AxesHandle(4),'Position',[dx,dy+dy/2,x,y]);
set(AxesHandle(3),'Position',[x+2.5*dx,dy+dy/2,x,y]);
set(AxesHandle(2),'Position',[2*x+4*dx,dy+dy/2,x,y]);
set(gcf,'Position',[100, 100, 1100, 700])



iter=iter_set;
figure5 = figure(5);
set(gcf,'color','w')

for i=1:NI-2
    if i==1
        ind = 2;
    elseif i>1 && i<4
        ind = 3;
%         col = (i-1)/3;
    elseif i>3 && i<6
        ind = 4;
%         col = (i-3)/4;
    elseif i>5 && i<8
        ind = 5;
%         col = (i-6)/5;
    elseif i>7 && i<10
        ind = 6;
%         col = (i-10)/6;
    elseif i == 10
        ind = 1;
%         col = (i-10)/6;
    end
    
    if i==2 || i==4 || i==6 || i==8
        flag_g = 1;
    elseif i==3 || i==5 || i==7 || i==9
        flag_g = 2;
    elseif i==10
        flag_g = 0;
    elseif i == 1
        flag_g = 3;
    else
        flag_g = 4;
    end
    
    refer_use =  bid_price_minus_ta{iter,1}{i,1}(t+1);
%     refer_use = 0;
   

    subplot(2,3,ind)

    h10 = plot([-t,t],[0,0],':','color',[1/2 1/2 1/2]);
    hold on
    
%     co = [0.0336;0.0380;0.6581;0.0493;0.8230;0.0627;0.9409;0.0772;1.0346;0];
    co = [-0.0307;-0.0338;-0.5384;-0.0439;-0.6491;-0.0556;-0.7197;-0.0694;-0.7625;0];
    
    if flag_g == 1
        t1 = plot(1,co(i,1),'s','color',[0 0 0],'LineWidth',1,'MarkerSize',10);
    elseif flag_g == 2
        t2 = plot(1,co(i,1),'x','color',[0 0 0],'LineWidth',1,'MarkerSize',10);
    elseif flag_g == 0
        t3 = plot(1,co(i,1),'o','color',[0 0 0],'LineWidth',1,'MarkerSize',10);
    elseif flag_g == 3
        t3 = plot(1,co(i,1),'d','color',[0 0 0],'LineWidth',1,'MarkerSize',10);
    end
    hold on

    if flag_g==1 % type A
        h1 = plot(-t:t,bid_price_minus_ta{iter,1}{i,1} - refer_use,'-','color',C(2,:),'LineWidth',1);
        bb = bid_price_minus_ta{iter,1}{i,1} - refer_use + 3;
        hold on
        h2 = plot(-t:t,ask_price_minus_ta{iter,1}{i,1} - refer_use,'-','color',C(1,:),'LineWidth',1);
        aa = ask_price_minus_ta{iter,1}{i,1} - refer_use + 3;
        hold on
%         plot(0,bid_price_minus_ta{iter,1}{i,1}(t+1) - refer_use,'*','color',C(2,:),'MarkerSize',10)
%         hold on
%         plot(0,ask_price_minus_ta{iter,1}{i,1}(t+1) - refer_use,'*','color',C(1,:),'MarkerSize',10)
%         hold on
%         plot(1,bid_price_minus_ta{iter,1}{i,1}(t+2) - refer_use,'*','color',C(2,:),'MarkerSize',10)
%         hold on
%         plot(1,ask_price_minus_ta{iter,1}{i,1}(t+2) - refer_use,'*','color',C(1,:),'MarkerSize',10)
    elseif flag_g==2 % type B
        h3 = plot(-t:t,bid_price_minus_ta{iter,1}{i,1} - refer_use,'--','color',C(2,:),'LineWidth',1);
        bb = bid_price_minus_ta{iter,1}{i,1} - refer_use + 3;
        hold on
        h4 = plot(-t:t,ask_price_minus_ta{iter,1}{i,1} - refer_use,'--','color',C(1,:),'LineWidth',1);
        aa = ask_price_minus_ta{iter,1}{i,1} - refer_use + 3;
        hold on
%         plot(-t:t,bid_price_minus_ta{iter,1}{i,1} - refer_use,'r.','LineWidth',1)
%         hold on
%         plot(-t:t,ask_price_minus_ta{iter,1}{i,1} - refer_use,'b.','LineWidth',1)
%         hold on
%         plot(0,bid_price_minus_ta{iter,1}{i,1}(t+1) - refer_use,'o','color',C(2,:),'MarkerSize',10)
%         hold on
%         plot(0,ask_price_minus_ta{iter,1}{i,1}(t+1) - refer_use,'o','color',C(1,:),'MarkerSize',10)
%         hold on
%         plot(1,bid_price_minus_ta{iter,1}{i,1}(t+2) - refer_use,'o','color',C(2,:),'MarkerSize',10)
%         hold on
%         plot(1,ask_price_minus_ta{iter,1}{i,1}(t+2) - refer_use,'o','color',C(1,:),'MarkerSize',10)
    elseif flag_g == 0% type Zero
        h5 = plot(-t:t,bid_price_minus_ta{iter,1}{i,1} - refer_use,'-.','color',C(2,:),'LineWidth',1);
        bb = bid_price_minus_ta{iter,1}{i,1} - refer_use + 3;
        hold on
        h6 = plot(-t:t,ask_price_minus_ta{iter,1}{i,1} - refer_use,'-.','color',C(1,:),'LineWidth',1);
        aa = ask_price_minus_ta{iter,1}{i,1} - refer_use + 3;
        hold on
    elseif flag_g == 3% type One
        h7 = plot(-t:t,bid_price_minus_ta{iter,1}{i,1} - refer_use,'-.','color',C(2,:),'LineWidth',1);
        bb = bid_price_minus_ta{iter,1}{i,1} - refer_use + 3;
        hold on
        h8 = plot(-t:t,ask_price_minus_ta{iter,1}{i,1} - refer_use,'-.','color',C(1,:),'LineWidth',1);
        aa = ask_price_minus_ta{iter,1}{i,1} - refer_use + 3;
        hold on
%         plot(0,bid_price_minus_ta{iter,1}{i,1}(t+1) - refer_use,'^','color',C(2,:),'MarkerSize',10)
%         hold on
%         plot(0,ask_price_minus_ta{iter,1}{i,1}(t+1) - refer_use,'^','color',C(1,:),'MarkerSize',10)
%         hold on
%         plot(1,bid_price_minus_ta{iter,1}{i,1}(t+2) - refer_use,'^','color',C(2,:),'MarkerSize',10)
%         hold on
%         plot(1,ask_price_minus_ta{iter,1}{i,1}(t+2) - refer_use,'^','color',C(1,:),'MarkerSize',10)
    end
    ylim([-5 3])
    xlim([-delta_t delta_t])
    ra_aa(i) = mean(aa(end-49:end)) - mean(aa(1:50));
    ra_bb(i) = mean(bb(end-49:end)) - mean(bb(1:50));
%     ylim([-2.5 5.2])
%     xlim([-.5 100])
%         set(gca,'xscale','log')
%         set(gca,'yscale','log')
%         hold on
%     hold on
%     plot([-t:1:-1,1:t],spread_minus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_minus_num{iter,1}(i,1),'b.')
%     legend('minus initiated market order','plus initiated market order')
    xlabel('\tau')
    if i == 1 || i == 3 || i== 5 || i == 7 || i == 9 || i == 10
        if ind==2
            ylabel('G_a(\tau|\Deltab_0=1),  G_b(\tau|\Deltab_0=1)')
            lh = legend([h7 h8],{'Bid (Type One)','Ask (Type One)'});
            legend boxoff
            set(lh, 'Position', [.197 .92-.28-.04 .5 .01])
        elseif ind==3
            ylabel('G_a(\tau|\Deltab_0=2),  G_b(\tau|\Deltab_0=2)')
            lh = legend([h1 h2 h3 h4],{'Bid (Type A)','Ask (Type A)','Bid (Type B)','Ask (Type B)'});
            legend boxoff
            set(lh, 'Position', [.53 .90-.28 .5 .01])
        elseif ind==4
            ylabel('G_a(\tau|\Deltab_0=3),  G_b(\tau|\Deltab_0=3)')
            lh = legend([h1 h2 h3 h4],{'Bid (Type A)','Ask (Type A)','Bid (Type B)','Ask (Type B)'});
            legend boxoff
            set(lh, 'Position', [-.135 .41-.28 .5 .01])
        elseif ind==5
            ylabel('G_a(\tau|\Deltab_0=4),  G_b(\tau|\Deltab_0=4)')
            lh = legend([h1 h2 h3 h4],{'Bid (Type A)','Ask (Type A)','Bid (Type B)','Ask (Type B)'});
            legend boxoff
            set(lh, 'Position', [.197 .41-.28 .5 .01])
        elseif ind==6
            ylabel('G_a(\tau|\Deltab_0=5),  G_b(\tau|\Deltab_0=5)')
            lh = legend([h1 h2 h3 h4],{'Bid (Type A)','Ask (Type A)','Bid (Type B)','Ask (Type B)'});
            legend boxoff
            set(lh, 'Position', [.53 .41-.28 .5 .01])
        elseif ind==1
            ylabel('G_a(\tau|\Deltab_0=0),  G_b(\tau|\Deltab_0=0)')
            lh = legend([h5 h6],{'Bid (Type Zero)','Ask (Type Zero)'});
            legend boxoff
            set(lh, 'Position', [-.128 .92-.28-.04 .5 .01])
        end
        if i==1
            text(-0.1,1.1,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
        elseif i==3
            text(-0.1,1.1,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top')
        elseif i==5
            text(-0.1,1.1,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top')
        elseif i==7
            text(-0.1,1.1,'(e)','Units', 'Normalized', 'VerticalAlignment', 'Top')
        elseif i==9
            text(-0.1,1.1,'(f)','Units', 'Normalized', 'VerticalAlignment', 'Top')
        elseif i==10
            text(-0.1,1.1,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
        end
    end
end
dx=0.05;
dy=0.05;
x=(1-5*dx)/2.9;
y=(1-5*dy)/1.9;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure5,'Type','axes');
set(AxesHandle(1),'Position',[dx,y+3.45*dy,x,y]);
set(AxesHandle(6),'Position',[x+2.5*dx,y+3.45*dy,x,y]);
set(AxesHandle(5),'Position',[2*x+4*dx,y+3.45*dy,x,y]);
set(AxesHandle(4),'Position',[dx,dy+dy/2,x,y]);
set(AxesHandle(3),'Position',[x+2.5*dx,dy+dy/2,x,y]);
set(AxesHandle(2),'Position',[2*x+4*dx,dy+dy/2,x,y]);
set(gcf,'Position',[100, 100, 1100, 700])



% % % % iter=iter_set;
% % % % f2=figure;
% % % % set(gcf,'color','w')
% % % % % set(f1,'Visible','off')
% % % % % set(f1, 'PaperUnits', 'inches');
% % % % % x_width=14 ;y_width=7;
% % % % % set(f1, 'PaperPosition', [0 0 x_width y_width]);
% % % % 
% % % % 
% % % % for i=1:NI
% % % %     if i>5 && i<10
% % % %         ind = i-4;
% % % %         flag = 1;
% % % %         flag_g = 1;
% % % %     elseif i>9
% % % %         ind = i-8;
% % % %         flag = 1;
% % % %         flag_g = 2;
% % % %     elseif i==1
% % % %         ind = 1;
% % % %         flag = 1;
% % % %         flag_g = 1;
% % % %     elseif i>1 && i<6
% % % %         flag = 0;
% % % %     end
% % % %     
% % % %     refer_use =  bid_price_minus_ta{iter,1}{i,1}(t+1);
% % % % %     refer_use = 0;
% % % %     
% % % %     if flag==1
% % % %         subplot(1,5,ind)
% % % %         
% % % %         if flag_g==1
% % % %             plot(-t:t,bid_price_minus_ta{iter,1}{i,1} - refer_use,'r-','LineWidth',1)
% % % %             hold on
% % % %             plot(-t:t,ask_price_minus_ta{iter,1}{i,1} - refer_use,'b-','LineWidth',1)
% % % %             hold on
% % % %             plot(0,bid_price_minus_ta{iter,1}{i,1}(t+1) - refer_use,'r*','MarkerSize',10)
% % % %             hold on
% % % %             plot(0,ask_price_minus_ta{iter,1}{i,1}(t+1) - refer_use,'b*','MarkerSize',10)
% % % %         elseif flag_g==2
% % % %             plot(-t:t,bid_price_minus_ta{iter,1}{i,1} - refer_use,'r.','LineWidth',1)
% % % %             hold on
% % % %             plot(-t:t,ask_price_minus_ta{iter,1}{i,1} - refer_use,'b.','LineWidth',1)
% % % %             hold on
% % % %             plot(0,bid_price_minus_ta{iter,1}{i,1}(t+1) - refer_use,'ro','MarkerSize',10)
% % % %             hold on
% % % %             plot(0,ask_price_minus_ta{iter,1}{i,1}(t+1) - refer_use,'bo','MarkerSize',10)
% % % %         end
% % % %         ylim([-8 3])
% % % % %         set(gca,'xscale','log')
% % % % %         set(gca,'yscale','log')
% % % % %         hold on
% % % %     %     hold on
% % % %     %     plot([-t:1:-1,1:t],spread_minus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_minus_num{iter,1}(i,1),'b.')
% % % %     %     legend('minus initiated market order','plus initiated market order')
% % % %         xlabel('lag')
% % % %         ylabel('Price')
% % % %     end
% % % % end
% % % 


% immediate & permanent impact conditioned on the ask price change
for i=1:NI
%     if i==1
%         ind = 1;
%     elseif i>1 && i<4
%         ind = 2;
% %         col = (i-1)/3;
%     elseif i>3 && i<6
%         ind = 3;
% %         col = (i-3)/4;
%     elseif i>5 && i<8
%         ind = 4;
% %         col = (i-6)/5;
%     elseif i>7 && i<10
%         ind = 5;
% %         col = (i-10)/6;
%     end
%     
%     if i==1 || i==2 || i==4 || i==6 || i==8
%         flag_g = 1;
%     elseif i==3 || i==5 || i==7 || i==9
%         flag_g = 2;
%     else
%         flag_g = 0;
%     end
    
    refer_use =  ask_price_plus_ta{iter,1}{i,1}(t+1);
    refer_use_mid =  (bid_price_plus_ta{iter,1}{i,1}(t+1) + ask_price_plus_ta{iter,1}{i,1}(t+1))/2;
%     refer_use = 0;
    
    bid_side = bid_price_plus_ta{iter,1}{i,1} - refer_use;
    ask_side = ask_price_plus_ta{iter,1}{i,1} - refer_use;
    spread = spread_plus_ta{iter,1}{i,1};
    
    mid = (bid_price_plus_ta{iter,1}{i,1} + ask_price_plus_ta{iter,1}{i,1})/2 - refer_use_mid;
    
    bid_vol = bid_vol_plus_ta{iter,1}{i,1};
    ask_vol = ask_vol_plus_ta{iter,1}{i,1};
    
    bid_gap = bid_gap_plus_ta{iter,1}{i,1};
    ask_gap = ask_gap_plus_ta{iter,1}{i,1};
    
%     shock_ind = 3;
%     bid_shock = bid_vol_price_plus_avg_ta{iter,1}{i,shock_ind};
%     ask_shock = ask_vol_price_plus_avg_ta{iter,1}{i,shock_ind};
    
    imm_bid_side(i,1) = bid_side(t+2,1) - bid_side(t+1,1);
    imm_ask_side(i,1) = ask_side(t+2,1) - ask_side(t+1,1);
    per_bid_side(i,1) = bid_side(we,1) - bid_side(t+2,1);
    per_ask_side(i,1) = ask_side(we,1) - ask_side(t+2,1);
    imm_spread_side(i,1) = spread(t+2,1) - spread(t+1,1);
    per_spread_side(i,1) = spread(we,1) - spread(t+2,1);
    imm_mid_side(i,1) = mid(t+2,1) - mid(t+1,1);
    per_mid_side(i,1) = mid(we,1) - mid(t+2,1);
    
    imm_bid_vol(i,1) = bid_vol(t+2,1) - bid_vol(t+1,1);
    imm_ask_vol(i,1) = ask_vol(t+2,1) - ask_vol(t+1,1);
    per_bid_vol(i,1) = bid_vol(we,1) - bid_vol(t+2,1);
    per_ask_vol(i,1) = ask_vol(we,1) - ask_vol(t+2,1);
    
    imm_bid_gap(i,1) = bid_gap(t+2,1) - bid_gap(t+1,1);
    imm_ask_gap(i,1) = ask_gap(t+2,1) - ask_gap(t+1,1);
    per_bid_gap(i,1) = bid_gap(we,1) - bid_gap(t+2,1);
    per_ask_gap(i,1) = ask_gap(we,1) - ask_gap(t+2,1);
    
%     imm_bid_shock(i,1) = bid_shock(t+2,1) - bid_shock(t+1,1);
%     imm_ask_shock(i,1) = ask_shock(t+2,1) - ask_shock(t+1,1);
%     per_bid_shock(i,1) = bid_shock(we,1) - bid_shock(t+2,1);
%     per_ask_shock(i,1) = ask_shock(we,1) - ask_shock(t+2,1);
    
    clear bid_side ask_side spread
    clear mid bid_vol ask_vol bid_gap ask_gap bid_shock ask_shock
    
    
    refer_use_minus =  bid_price_minus_ta{iter,1}{i,1}(t+1);
    refer_use_mid_minus =  (bid_price_minus_ta{iter,1}{i,1}(t+1) + ask_price_minus_ta{iter,1}{i,1}(t+1))/2;
%     refer_use = 0;
    
    bid_side_minus = bid_price_minus_ta{iter,1}{i,1} - refer_use_minus;
    ask_side_minus = ask_price_minus_ta{iter,1}{i,1} - refer_use_minus;
    spread_minus = spread_minus_ta{iter,1}{i,1};
    
    mid_minus = (bid_price_minus_ta{iter,1}{i,1} + ask_price_minus_ta{iter,1}{i,1})/2 - refer_use_mid_minus;
    
    bid_vol_minus = bid_vol_minus_ta{iter,1}{i,1};
    ask_vol_minus = ask_vol_minus_ta{iter,1}{i,1};
    
    bid_gap_minus = bid_gap_minus_ta{iter,1}{i,1};
    ask_gap_minus = ask_gap_minus_ta{iter,1}{i,1};
    
%     shock_ind = 3;
%     bid_shock_minus = bid_vol_price_minus_avg_ta{iter,1}{i,shock_ind};
%     ask_shock_minus = ask_vol_price_minus_avg_ta{iter,1}{i,shock_ind};
    
    imm_bid_side_minus(i,1) = bid_side_minus(t+2,1) - bid_side_minus(t+1,1);
    imm_ask_side_minus(i,1) = ask_side_minus(t+2,1) - ask_side_minus(t+1,1);
    per_bid_side_minus(i,1) = bid_side_minus(we,1) - bid_side_minus(t+2,1);
    per_ask_side_minus(i,1) = ask_side_minus(we,1) - ask_side_minus(t+2,1);
    imm_spread_side_minus(i,1) = spread_minus(t+2,1) - spread_minus(t+1,1);
    per_spread_side_minus(i,1) = spread_minus(we,1) - spread_minus(t+2,1);
    imm_mid_side_minus(i,1) = mid_minus(t+2,1) - mid_minus(t+1,1);
    per_mid_side_minus(i,1) = mid_minus(we,1) - mid_minus(t+2,1);
    
    imm_bid_vol_minus(i,1) = bid_vol_minus(t+2,1) - bid_vol_minus(t+1,1);
    imm_ask_vol_minus(i,1) = ask_vol_minus(t+2,1) - ask_vol_minus(t+1,1);
    per_bid_vol_minus(i,1) = bid_vol_minus(we,1) - bid_vol_minus(t+2,1);
    per_ask_vol_minus(i,1) = ask_vol_minus(we,1) - ask_vol_minus(t+2,1);
    
    imm_bid_gap_minus(i,1) = bid_gap_minus(t+2,1) - bid_gap_minus(t+1,1);
    imm_ask_gap_minus(i,1) = ask_gap_minus(t+2,1) - ask_gap_minus(t+1,1);
    per_bid_gap_minus(i,1) = bid_gap_minus(we,1) - bid_gap_minus(t+2,1);
    per_ask_gap_minus(i,1) = ask_gap_minus(we,1) - ask_gap_minus(t+2,1);
    
%     imm_bid_shock_minus(i,1) = bid_shock_minus(t+2,1) - bid_shock_minus(t+1,1);
%     imm_ask_shock_minus(i,1) = ask_shock_minus(t+2,1) - ask_shock_minus(t+1,1);
%     per_bid_shock_minus(i,1) = bid_shock_minus(we,1) - bid_shock_minus(t+2,1);
%     per_ask_shock_minus(i,1) = ask_shock_minus(we,1) - ask_shock_minus(t+2,1);
    
    clear bid_side_minus ask_side_minus spread_minus
    clear mid_minus bid_vol_minus ask_vol_minus bid_gap_minus ask_gap_minus bid_shock_minus ask_shock_minus

%     subplot(1,5,ind)

%     if flag_g==1
%         subplot(1,2,1)
%         scatter(imm_bid_side(i,1),per_bid_side(i,1),[],[1 col col])
%         hold on
%         subplot(1,2,2)
%         scatter(imm_ask_side(i,1),per_ask_side(i,1),[],[col col 1])
%         hold on
%     elseif flag_g==2
%         subplot(1,2,1)
%         scatter(imm_bid_side(i,1),per_bid_side(i,1),[],[1 col col])
%         hold on
%         subplot(1,2,2)
%         scatter(imm_ask_side(i,1),per_ask_side(i,1),[],[col col 1])
%         hold on
%     end
%     xlabel('lag')
%     ylabel('Price')
end

% immediate
imm_bid(:,1) = [imm_bid_side(1); imm_bid_side(2); imm_bid_side(4); imm_bid_side(6); imm_bid_side(8)];
imm_ask(:,1) = [imm_ask_side(1); imm_ask_side(2); imm_ask_side(4); imm_ask_side(6); imm_ask_side(8)];
imm_bid(:,2) = [imm_bid_side(1); imm_bid_side(3); mean(imm_bid_side(5)); mean(imm_bid_side(7)); mean(imm_bid_side(9))];
imm_ask(:,2) = [imm_ask_side(1); imm_ask_side(3); mean(imm_ask_side(5)); mean(imm_ask_side(7)); mean(imm_ask_side(9))];
imm_spread(:,1) = [imm_spread_side(1); imm_spread_side(2); imm_spread_side(4); imm_spread_side(6); imm_spread_side(8)];
imm_spread(:,2) = [imm_spread_side(1); imm_spread_side(3); mean(imm_spread_side(5)); mean(imm_spread_side(7)); mean(imm_spread_side(9))];
imm_mid(:,1) = [imm_mid_side(1); imm_mid_side(2); imm_mid_side(4); imm_mid_side(6); imm_mid_side(8)];
imm_mid(:,2) = [imm_mid_side(1); imm_mid_side(3); mean(imm_mid_side(5)); mean(imm_mid_side(7)); mean(imm_mid_side(9))];

imm_bid_vol_result(:,1) = [imm_bid_vol(1); imm_bid_vol(2); imm_bid_vol(4); imm_bid_vol(6); imm_bid_vol(8)];
imm_ask_vol_result(:,1) = [imm_ask_vol(1); imm_ask_vol(2); imm_ask_vol(4); imm_ask_vol(6); imm_ask_vol(8)];
imm_bid_vol_result(:,2) = [imm_bid_vol(1); imm_bid_vol(3); mean(imm_bid_vol(5)); mean(imm_bid_vol(7)); mean(imm_bid_vol(9))];
imm_ask_vol_result(:,2) = [imm_ask_vol(1); imm_ask_vol(3); mean(imm_ask_vol(5)); mean(imm_ask_vol(7)); mean(imm_ask_vol(9))];
imm_bid_gap_result(:,1) = [imm_bid_gap(1); imm_bid_gap(2); imm_bid_gap(4); imm_bid_gap(6); imm_bid_gap(8)];
imm_ask_gap_result(:,1) = [imm_ask_gap(1); imm_ask_gap(2); imm_ask_gap(4); imm_ask_gap(6); imm_ask_gap(8)];
imm_bid_gap_result(:,2) = [imm_bid_gap(1); imm_bid_gap(3); mean(imm_bid_gap(5)); mean(imm_bid_gap(7)); mean(imm_bid_gap(9))];
imm_ask_gap_result(:,2) = [imm_ask_gap(1); imm_ask_gap(3); mean(imm_ask_gap(5)); mean(imm_ask_gap(7)); mean(imm_ask_gap(9))];
% imm_bid_shock_result(:,1) = [imm_bid_shock(1); imm_bid_shock(2); imm_bid_shock(4); imm_bid_shock(6); imm_bid_shock(8)];
% imm_ask_shock_result(:,1) = [imm_ask_shock(1); imm_ask_shock(2); imm_ask_shock(4); imm_ask_shock(6); imm_ask_shock(8)];
% imm_bid_shock_result(:,2) = [imm_bid_shock(1); imm_bid_shock(3); mean(imm_bid_shock(5)); mean(imm_bid_shock(7)); mean(imm_bid_shock(9))];
% imm_ask_shock_result(:,2) = [imm_ask_shock(1); imm_ask_shock(3); mean(imm_ask_shock(5)); mean(imm_ask_shock(7)); mean(imm_ask_shock(9))];


% reversion
per_bid(:,1) = [per_bid_side(1); per_bid_side(2); per_bid_side(4); per_bid_side(6); per_bid_side(8)];
per_ask(:,1) = [per_ask_side(1); per_ask_side(2); per_ask_side(4); per_ask_side(6); per_ask_side(8)];
per_bid(:,2) = [per_bid_side(1); per_bid_side(3); mean(per_bid_side(5)); mean(per_bid_side(7)); mean(per_bid_side(9))];
per_ask(:,2) = [per_ask_side(1); per_ask_side(3); mean(per_ask_side(5)); mean(per_ask_side(7)); mean(per_ask_side(9))];
per_spread(:,1) = [per_spread_side(1); per_spread_side(2); per_spread_side(4); per_spread_side(6); per_spread_side(8)];
per_spread(:,2) = [per_spread_side(1); per_spread_side(3); mean(per_spread_side(5)); mean(per_spread_side(7)); mean(per_spread_side(9))];
per_mid(:,1) = [per_mid_side(1); per_mid_side(2); per_mid_side(4); per_mid_side(6); per_mid_side(8)];
per_mid(:,2) = [per_mid_side(1); per_mid_side(3); mean(per_mid_side(5)); mean(per_mid_side(7)); mean(per_mid_side(9))];

per_bid_vol_result(:,1) = [per_bid_vol(1); per_bid_vol(2); per_bid_vol(4); per_bid_vol(6); per_bid_vol(8)];
per_ask_vol_result(:,1) = [per_ask_vol(1); per_ask_vol(2); per_ask_vol(4); per_ask_vol(6); per_ask_vol(8)];
per_bid_vol_result(:,2) = [per_bid_vol(1); per_bid_vol(3); mean(per_bid_vol(5)); mean(per_bid_vol(7)); mean(per_bid_vol(9))];
per_ask_vol_result(:,2) = [per_ask_vol(1); per_ask_vol(3); mean(per_ask_vol(5)); mean(per_ask_vol(7)); mean(per_ask_vol(9))];
per_bid_gap_result(:,1) = [per_bid_gap(1); per_bid_gap(2); per_bid_gap(4); per_bid_gap(6); per_bid_gap(8)];
per_ask_gap_result(:,1) = [per_ask_gap(1); per_ask_gap(2); per_ask_gap(4); per_ask_gap(6); per_ask_gap(8)];
per_bid_gap_result(:,2) = [per_bid_gap(1); per_bid_gap(3); mean(per_bid_gap(5)); mean(per_bid_gap(7)); mean(per_bid_gap(9))];
per_ask_gap_result(:,2) = [per_ask_gap(1); per_ask_gap(3); mean(per_ask_gap(5)); mean(per_ask_gap(7)); mean(per_ask_gap(9))];
% per_bid_shock_result(:,1) = [per_bid_shock(1); per_bid_shock(2); per_bid_shock(4); per_bid_shock(6); per_bid_shock(8)];
% per_ask_shock_result(:,1) = [per_ask_shock(1); per_ask_shock(2); per_ask_shock(4); per_ask_shock(6); per_ask_shock(8)];
% per_bid_shock_result(:,2) = [per_bid_shock(1); per_bid_shock(3); mean(per_bid_shock(5)); mean(per_bid_shock(7)); mean(per_bid_shock(9))];
% per_ask_shock_result(:,2) = [per_ask_shock(1); per_ask_shock(3); mean(per_ask_shock(5)); mean(per_ask_shock(7)); mean(per_ask_shock(9))];




imm_bid_minus(:,1) = [imm_bid_side_minus(1); imm_bid_side_minus(2); imm_bid_side_minus(4); imm_bid_side_minus(6); imm_bid_side_minus(8)];
imm_ask_minus(:,1) = [imm_ask_side_minus(1); imm_ask_side_minus(2); imm_ask_side_minus(4); imm_ask_side_minus(6); imm_ask_side_minus(8)];
imm_bid_minus(:,2) = [imm_bid_side_minus(1); imm_bid_side_minus(3); mean(imm_bid_side_minus(5)); mean(imm_bid_side_minus(7)); mean(imm_bid_side_minus(9))];
imm_ask_minus(:,2) = [imm_ask_side_minus(1); imm_ask_side_minus(3); mean(imm_ask_side_minus(5)); mean(imm_ask_side_minus(7)); mean(imm_ask_side_minus(9))];
imm_spread_minus(:,1) = [imm_spread_side_minus(1); imm_spread_side_minus(2); imm_spread_side_minus(4); imm_spread_side_minus(6); imm_spread_side_minus(8)];
imm_spread_minus(:,2) = [imm_spread_side_minus(1); imm_spread_side_minus(3); mean(imm_spread_side_minus(5)); mean(imm_spread_side_minus(7)); mean(imm_spread_side_minus(9))];
imm_mid_minus(:,1) = [imm_mid_side_minus(1); imm_mid_side_minus(2); imm_mid_side_minus(4); imm_mid_side_minus(6); imm_mid_side_minus(8)];
imm_mid_minus(:,2) = [imm_mid_side_minus(1); imm_mid_side_minus(3); mean(imm_mid_side_minus(5)); mean(imm_mid_side_minus(7)); mean(imm_mid_side_minus(9))];
%             set(lh, 'Position', [.2 .95 .5 .01])
imm_bid_vol_minus_result(:,1) = [imm_bid_vol_minus(1); imm_bid_vol_minus(2); imm_bid_vol_minus(4); imm_bid_vol_minus(6); imm_bid_vol_minus(8)];
imm_ask_vol_minus_result(:,1) = [imm_ask_vol_minus(1); imm_ask_vol_minus(2); imm_ask_vol_minus(4); imm_ask_vol_minus(6); imm_ask_vol_minus(8)];
imm_bid_vol_minus_result(:,2) = [imm_bid_vol_minus(1); imm_bid_vol_minus(3); mean(imm_bid_vol_minus(5)); mean(imm_bid_vol_minus(7)); mean(imm_bid_vol_minus(9))];
imm_ask_vol_minus_result(:,2) = [imm_ask_vol_minus(1); imm_ask_vol_minus(3); mean(imm_ask_vol_minus(5)); mean(imm_ask_vol_minus(7)); mean(imm_ask_vol_minus(9))];
imm_bid_gap_minus_result(:,1) = [imm_bid_gap_minus(1); imm_bid_gap_minus(2); imm_bid_gap_minus(4); imm_bid_gap_minus(6); imm_bid_gap_minus(8)];
imm_ask_gap_minus_result(:,1) = [imm_ask_gap_minus(1); imm_ask_gap_minus(2); imm_ask_gap_minus(4); imm_ask_gap_minus(6); imm_ask_gap_minus(8)];
imm_bid_gap_minus_result(:,2) = [imm_bid_gap_minus(1); imm_bid_gap_minus(3); mean(imm_bid_gap_minus(5)); mean(imm_bid_gap_minus(7)); mean(imm_bid_gap_minus(9))];
imm_ask_gap_minus_result(:,2) = [imm_ask_gap_minus(1); imm_ask_gap_minus(3); mean(imm_ask_gap_minus(5)); mean(imm_ask_gap_minus(7)); mean(imm_ask_gap_minus(9))];
% imm_bid_shock_minus_result(:,1) = [imm_bid_shock_minus(1); imm_bid_shock_minus(2); imm_bid_shock_minus(4); imm_bid_shock_minus(6); imm_bid_shock_minus(8)];
% imm_ask_shock_minus_result(:,1) = [imm_ask_shock_minus(1); imm_ask_shock_minus(2); imm_ask_shock_minus(4); imm_ask_shock_minus(6); imm_ask_shock_minus(8)];
% imm_bid_shock_minus_result(:,2) = [imm_bid_shock_minus(1); imm_bid_shock_minus(3); mean(imm_bid_shock_minus(5)); mean(imm_bid_shock_minus(7)); mean(imm_bid_shock_minus(9))];
% imm_ask_shock_minus_result(:,2) = [imm_ask_shock_minus(1); imm_ask_shock_minus(3); mean(imm_ask_shock_minus(5)); mean(imm_ask_shock_minus(7)); mean(imm_ask_shock_minus(9))];


% reversion
per_bid_minus(:,1) = [per_bid_side_minus(1); per_bid_side_minus(2); per_bid_side_minus(4); per_bid_side_minus(6); per_bid_side_minus(8)];
per_ask_minus(:,1) = [per_ask_side_minus(1); per_ask_side_minus(2); per_ask_side_minus(4); per_ask_side_minus(6); per_ask_side_minus(8)];
per_bid_minus(:,2) = [per_bid_side_minus(1); per_bid_side_minus(3); mean(per_bid_side_minus(5)); mean(per_bid_side_minus(7)); mean(per_bid_side_minus(9))];
per_ask_minus(:,2) = [per_ask_side_minus(1); per_ask_side_minus(3); mean(per_ask_side_minus(5)); mean(per_ask_side_minus(7)); mean(per_ask_side_minus(9))];
per_spread_minus(:,1) = [per_spread_side_minus(1); per_spread_side_minus(2); per_spread_side_minus(4); per_spread_side_minus(6); per_spread_side_minus(8)];
per_spread_minus(:,2) = [per_spread_side_minus(1); per_spread_side_minus(3); mean(per_spread_side_minus(5)); mean(per_spread_side_minus(7)); mean(per_spread_side_minus(9))];
per_mid_minus(:,1) = [per_mid_side_minus(1); per_mid_side_minus(2); per_mid_side_minus(4); per_mid_side_minus(6); per_mid_side_minus(8)];
per_mid_minus(:,2) = [per_mid_side_minus(1); per_mid_side_minus(3); mean(per_mid_side_minus(5)); mean(per_mid_side_minus(7)); mean(per_mid_side_minus(9))];

per_bid_vol_minus_result(:,1) = [per_bid_vol_minus(1); per_bid_vol_minus(2); per_bid_vol_minus(4); per_bid_vol_minus(6); per_bid_vol_minus(8)];
per_ask_vol_minus_result(:,1) = [per_ask_vol_minus(1); per_ask_vol_minus(2); per_ask_vol_minus(4); per_ask_vol_minus(6); per_ask_vol_minus(8)];
per_bid_vol_minus_result(:,2) = [per_bid_vol_minus(1); per_bid_vol_minus(3); mean(per_bid_vol_minus(5)); mean(per_bid_vol_minus(7)); mean(per_bid_vol_minus(9))];
per_ask_vol_minus_result(:,2) = [per_ask_vol_minus(1); per_ask_vol_minus(3); mean(per_ask_vol_minus(5)); mean(per_ask_vol_minus(7)); mean(per_ask_vol_minus(9))];
per_bid_gap_minus_result(:,1) = [per_bid_gap_minus(1); per_bid_gap_minus(2); per_bid_gap_minus(4); per_bid_gap_minus(6); per_bid_gap_minus(8)];
per_ask_gap_minus_result(:,1) = [per_ask_gap_minus(1); per_ask_gap_minus(2); per_ask_gap_minus(4); per_ask_gap_minus(6); per_ask_gap_minus(8)];
per_bid_gap_minus_result(:,2) = [per_bid_gap_minus(1); per_bid_gap_minus(3); mean(per_bid_gap_minus(5)); mean(per_bid_gap_minus(7)); mean(per_bid_gap_minus(9))];
per_ask_gap_minus_result(:,2) = [per_ask_gap_minus(1); per_ask_gap_minus(3); mean(per_ask_gap_minus(5)); mean(per_ask_gap_minus(7)); mean(per_ask_gap_minus(9))];
% per_bid_shock_minus_result(:,1) = [per_bid_shock_minus(1); per_bid_shock_minus(2); per_bid_shock_minus(4); per_bid_shock_minus(6); per_bid_shock_minus(8)];
% per_ask_shock_minus_result(:,1) = [per_ask_shock_minus(1); per_ask_shock_minus(2); per_ask_shock_minus(4); per_ask_shock_minus(6); per_ask_shock_minus(8)];
% per_bid_shock_minus_result(:,2) = [per_bid_shock_minus(1); per_bid_shock_minus(3); mean(per_bid_shock_minus(5)); mean(per_bid_shock_minus(7)); mean(per_bid_shock_minus(9))];
% per_ask_shock_minus_result(:,2) = [per_ask_shock_minus(1); per_ask_shock_minus(3); mean(per_ask_shock_minus(5)); mean(per_ask_shock_minus(7)); mean(per_ask_shock_minus(9))];


% ni = 15
% imm_bid(:,1) = [imm_bid_side(1); imm_bid_side(2); imm_bid_side(4); imm_bid_side(7); imm_bid_side(11)];
% imm_ask(:,1) = [imm_ask_side(1); imm_ask_side(2); imm_ask_side(4); imm_ask_side(7); imm_ask_side(11)];
% imm_bid(:,2) = [imm_bid_side(1); imm_bid_side(3); mean(imm_bid_side(5:6)); mean(imm_bid_side(8:10)); mean(imm_bid_side(12:14))];
% imm_ask(:,2) = [imm_ask_side(1); imm_ask_side(3); mean(imm_ask_side(5:6)); mean(imm_ask_side(8:10)); mean(imm_ask_side(12:14))];
% imm_spread(:,1) = [imm_spread_side(1); imm_spread_side(2); imm_spread_side(4); imm_spread_side(7); imm_spread_side(11)];
% imm_spread(:,2) = [imm_spread_side(1); imm_spread_side(3); mean(imm_spread_side(5:6)); mean(imm_spread_side(8:10)); mean(imm_spread_side(12:14))];
% imm_mid(:,1) = [imm_mid_side(1); imm_mid_side(2); imm_mid_side(4); imm_mid_side(7); imm_mid_side(11)];
% imm_mid(:,2) = [imm_mid_side(1); imm_mid_side(3); mean(imm_mid_side(5:6)); mean(imm_mid_side(8:10)); mean(imm_mid_side(12:14))];
% 
% imm_bid_vol_result(:,1) = [imm_bid_vol(1); imm_bid_vol(2); imm_bid_vol(4); imm_bid_vol(7); imm_bid_vol(11)];
% imm_ask_vol_result(:,1) = [imm_ask_vol(1); imm_ask_vol(2); imm_ask_vol(4); imm_ask_vol(7); imm_ask_vol(11)];
% imm_bid_vol_result(:,2) = [imm_bid_vol(1); imm_bid_vol(3); mean(imm_bid_vol(5:6)); mean(imm_bid_vol(8:10)); mean(imm_bid_vol(12:14))];
% imm_ask_vol_result(:,2) = [imm_ask_vol(1); imm_ask_vol(3); mean(imm_ask_vol(5:6)); mean(imm_ask_vol(8:10)); mean(imm_ask_vol(12:14))];
% imm_bid_gap_result(:,1) = [imm_bid_gap(1); imm_bid_gap(2); imm_bid_gap(4); imm_bid_gap(7); imm_bid_gap(11)];
% imm_ask_gap_result(:,1) = [imm_ask_gap(1); imm_ask_gap(2); imm_ask_gap(4); imm_ask_gap(7); imm_ask_gap(11)];
% imm_bid_gap_result(:,2) = [imm_bid_gap(1); imm_bid_gap(3); mean(imm_bid_gap(5:6)); mean(imm_bid_gap(8:10)); mean(imm_bid_gap(12:14))];
% imm_ask_gap_result(:,2) = [imm_ask_gap(1); imm_ask_gap(3); mean(imm_ask_gap(5:6)); mean(imm_ask_gap(8:10)); mean(imm_ask_gap(12:14))];
% imm_bid_shock_result(:,1) = [imm_bid_shock(1); imm_bid_shock(2); imm_bid_shock(4); imm_bid_shock(7); imm_bid_shock(11)];
% imm_ask_shock_result(:,1) = [imm_ask_shock(1); imm_ask_shock(2); imm_ask_shock(4); imm_ask_shock(7); imm_ask_shock(11)];
% imm_bid_shock_result(:,2) = [imm_bid_shock(1); imm_bid_shock(3); mean(imm_bid_shock(5:6)); mean(imm_bid_shock(8:10)); mean(imm_bid_shock(12:14))];
% imm_ask_shock_result(:,2) = [imm_ask_shock(1); imm_ask_shock(3); mean(imm_ask_shock(5:6)); mean(imm_ask_shock(8:10)); mean(imm_ask_shock(12:14))];
% 
% 
% % reversion
% per_bid(:,1) = [per_bid_side(1); per_bid_side(2); per_bid_side(4); per_bid_side(7); per_bid_side(11)];
% per_ask(:,1) = [per_ask_side(1); per_ask_side(2); per_ask_side(4); per_ask_side(7); per_ask_side(11)];
% per_bid(:,2) = [per_bid_side(1); per_bid_side(3); mean(per_bid_side(5:6)); mean(per_bid_side(8:10)); mean(per_bid_side(12:14))];
% per_ask(:,2) = [per_ask_side(1); per_ask_side(3); mean(per_ask_side(5:6)); mean(per_ask_side(8:10)); mean(per_ask_side(12:14))];
% per_spread(:,1) = [per_spread_side(1); per_spread_side(2); per_spread_side(4); per_spread_side(7); per_spread_side(11)];
% per_spread(:,2) = [per_spread_side(1); per_spread_side(3); mean(per_spread_side(5:6)); mean(per_spread_side(8:10)); mean(per_spread_side(12:14))];
% per_mid(:,1) = [per_mid_side(1); per_mid_side(2); per_mid_side(4); per_mid_side(7); per_mid_side(11)];
% per_mid(:,2) = [per_mid_side(1); per_mid_side(3); mean(per_mid_side(5:6)); mean(per_mid_side(8:10)); mean(per_mid_side(12:14))];
% 
% per_bid_vol_result(:,1) = [per_bid_vol(1); per_bid_vol(2); per_bid_vol(4); per_bid_vol(7); per_bid_vol(11)];
% per_ask_vol_result(:,1) = [per_ask_vol(1); per_ask_vol(2); per_ask_vol(4); per_ask_vol(7); per_ask_vol(11)];
% per_bid_vol_result(:,2) = [per_bid_vol(1); per_bid_vol(3); mean(per_bid_vol(5:6)); mean(per_bid_vol(8:10)); mean(per_bid_vol(12:14))];
% per_ask_vol_result(:,2) = [per_ask_vol(1); per_ask_vol(3); mean(per_ask_vol(5:6)); mean(per_ask_vol(8:10)); mean(per_ask_vol(12:14))];
% per_bid_gap_result(:,1) = [per_bid_gap(1); per_bid_gap(2); per_bid_gap(4); per_bid_gap(7); per_bid_gap(11)];
% per_ask_gap_result(:,1) = [per_ask_gap(1); per_ask_gap(2); per_ask_gap(4); per_ask_gap(7); per_ask_gap(11)];
% per_bid_gap_result(:,2) = [per_bid_gap(1); per_bid_gap(3); mean(per_bid_gap(5:6)); mean(per_bid_gap(8:10)); mean(per_bid_gap(12:14))];
% per_ask_gap_result(:,2) = [per_ask_gap(1); per_ask_gap(3); mean(per_ask_gap(5:6)); mean(per_ask_gap(8:10)); mean(per_ask_gap(12:14))];
% per_bid_shock_result(:,1) = [per_bid_shock(1); per_bid_shock(2); per_bid_shock(4); per_bid_shock(7); per_bid_shock(11)];
% per_ask_shock_result(:,1) = [per_ask_shock(1); per_ask_shock(2); per_ask_shock(4); per_ask_shock(7); per_ask_shock(11)];
% per_bid_shock_result(:,2) = [per_bid_shock(1); per_bid_shock(3); mean(per_bid_shock(5:6)); mean(per_bid_shock(8:10)); mean(per_bid_shock(12:14))];
% per_ask_shock_result(:,2) = [per_ask_shock(1); per_ask_shock(3); mean(per_ask_shock(5:6)); mean(per_ask_shock(8:10)); mean(per_ask_shock(12:14))];

% figure100 = figure('Position',[1 1 1000 1000],'Color',[1 1 1]);


% C2 = [204,37,41;57,106,177]/255;

figure2 = figure(2);
set(gcf,'color','w')
subplot(2,2,1)
hold on
% plot(imm_ask(:,1),imm_bid(:,1),'*-','color',C2(2,:),'LineWidth',lw);
h1 = plot(imm_ask(2:end,1),imm_bid(2:end,1),'*-','color',C2(2,:),'LineWidth',lw);
hold on
plot(imm_ask(1,1),imm_bid(1,1),'d','color',C2(2,:),'LineWidth',lw);
hold on
plot(imm_ask(1:2,1),imm_bid(1:2,1),'-','color',C2(2,:),'LineWidth',lw);
hold on
plot(imm_ask(1:2,2),imm_bid(1:2,2),'--','color',C2(2,:),'LineWidth',lw)
hold on
% plot(imm_ask(2:end,2),imm_bid(2:end,2),'o--','color',C2(2,:),'LineWidth',lw);
h2 = plot(imm_ask(2:end,2),imm_bid(2:end,2),'o--','color',C2(2,:),'LineWidth',lw);
hold on
plot(imm_ask(2:end,1),per_bid(2:end,1),'*-','color',C2(1,:),'LineWidth',lw);
hold on
plot(imm_ask(1:2,1),per_bid(1:2,1),'-','color',C2(1,:),'LineWidth',lw);
hold on
plot(imm_ask(1,1),per_bid(1,1),'d','color',C2(1,:),'LineWidth',lw);
hold on
plot(imm_ask(1:2,2),per_bid(1:2,2),'--','color',C2(1,:),'LineWidth',lw)
hold on
plot(imm_ask(2:end,2),per_bid(2:end,2),'o--','color',C2(1,:),'LineWidth',lw)
hold on
plot(0:1,[imm_bid_side(10);imm_bid(1,1)],'-.','color',C2(2,:),'LineWidth',lw)
hold on
plot(0:1,[per_bid_side(10);per_bid(1,1)],'-.','color',C2(1,:),'LineWidth',lw)
hold on
plot(0,[imm_bid_side(10)],'^-.','color',C2(2,:))
hold on
plot(0,[per_bid_side(10)],'^-.','color',C2(1,:))

% hold on
% plot(imm_bid_minus(:,1),imm_bid_minus(:,1),'k*-','LineWidth',lw)
% hold on
% plot(imm_bid_minus(1:2,2),imm_bid_minus(1:2,2),'k--','LineWidth',lw)
% hold on
% plot(imm_bid_minus(2:end,2),imm_bid_minus(2:end,2),'ko--','LineWidth',lw)
hold on
h3 = plot(imm_bid_minus(2:end,1),per_bid_minus(2:end,1),'*-','color',C2(1,:),'LineWidth',lw);
hold on
plot(imm_bid_minus(1:2,1),per_bid_minus(1:2,1),'-','color',C2(1,:),'LineWidth',lw);
hold on
plot(imm_bid_minus(1,1),per_bid_minus(1,1),'d','color',C2(1,:),'LineWidth',lw);
hold on
plot(imm_bid_minus(1:2,2),per_bid_minus(1:2,2),'--','color',C2(1,:),'LineWidth',lw)
hold on
h4 = plot(imm_bid_minus(2:end,2),per_bid_minus(2:end,2),'o--','color',C2(1,:),'LineWidth',lw);
hold on
% plot([0 -1],[imm_bid_side_minus(10);imm_bid_minus(1,1)],'k-.','LineWidth',lw)
hold on
plot([0 -1],[per_bid_side_minus(10);per_bid_minus(1,1)],'-.','color',C2(1,:),'LineWidth',lw)
hold on
% plot(0,[imm_bid_side_minus(10)],'k^-.')
hold on
plot(0,[per_bid_side_minus(10)],'^-.','color',C2(1,:))
lh = legend([h1 h2 h3 h4],'BID immediate change (type A)','BID immediate change (type B)','BID reversion (type A)','BID reversion (type B)');
legend boxoff
if lind == 1
set(lh, 'Position', [.01 .885 .5 .01])

xlabel('\DeltaA_0')
ylabel('E[R_b(\tau)|\DeltaA_0],   E[\Deltab_0|\DeltaA_0]')
text(-0.1,1.1,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlim([-5.5 5.5])
grid off
grid off
box on
end

subplot(2,2,2)
% plot(imm_ask(:,1),imm_ask(:,1),'k*-','LineWidth',lw)
% hold on
% plot(imm_ask(1:2,2),imm_ask(1:2,2),'k--','LineWidth',lw)
% hold on
% plot(imm_ask(2:end,2),imm_ask(2:end,2),'ko--','LineWidth',lw)
hold on
h3 = plot(imm_ask(2:end,1),per_ask(2:end,1),'*-','color',C2(1,:),'LineWidth',lw);
hold on
plot(imm_ask(1:2,1),per_ask(1:2,1),'-','color',C2(1,:),'LineWidth',lw);
hold on
plot(imm_ask(1,1),per_ask(1,1),'d','color',C2(1,:),'LineWidth',lw);
hold on
plot(imm_ask(1:2,2),per_ask(1:2,2),'--','color',C2(1,:),'LineWidth',lw)
hold on
h4 = plot(imm_ask(2:end,2),per_ask(2:end,2),'o--','color',C2(1,:),'LineWidth',lw);
hold on
% plot(0:1,[imm_ask_side(10);imm_ask(1,1)],'k-.','LineWidth',lw)
hold on
plot(0:1,[per_ask_side(10);per_ask(1,1)],'-.','color',C2(1,:),'LineWidth',lw)
hold on
% plot(0,[imm_ask_side(10)],'k^-.')
hold on
plot(0,[per_ask_side(10)],'^-.','color',C2(1,:))

hold on
h1 = plot(imm_bid_minus(2:end,1),imm_ask_minus(2:end,1),'*-','color',C2(2,:),'LineWidth',lw);
hold on
plot(imm_bid_minus(1:2,1),imm_ask_minus(1:2,1),'-','color',C2(2,:),'LineWidth',lw);
hold on
plot(imm_bid_minus(1,1),imm_ask_minus(1,1),'d','color',C2(2,:),'LineWidth',lw);
hold on
plot(imm_bid_minus(1:2,2),imm_ask_minus(1:2,2),'--','color',C2(2,:),'LineWidth',lw)
hold on
h2 = plot(imm_bid_minus(2:end,2),imm_ask_minus(2:end,2),'o--','color',C2(2,:),'LineWidth',lw);
hold on
plot(imm_bid_minus(2:end,1),per_ask_minus(2:end,1),'*-','color',C2(1,:),'LineWidth',lw)
hold on
plot(imm_bid_minus(1:2,1),per_ask_minus(1:2,1),'-','color',C2(1,:),'LineWidth',lw)
hold on
plot(imm_bid_minus(1,1),per_ask_minus(1,1),'d','color',C2(1,:),'LineWidth',lw)
hold on
plot(imm_bid_minus(1:2,2),per_ask_minus(1:2,2),'--','color',C2(1,:),'LineWidth',lw)
hold on
plot(imm_bid_minus(2:end,2),per_ask_minus(2:end,2),'o--','color',C2(1,:),'LineWidth',lw)
hold on
plot([0 -1],[imm_ask_side_minus(10);imm_ask_minus(1,1)],'-.','color',C2(2,:),'LineWidth',lw)
hold on
plot([0 -1],[per_ask_side_minus(10);per_ask_minus(1,1)],'-.','color',C2(1,:),'LineWidth',lw)
hold on
plot(0,[imm_ask_side_minus(10)],'^-.','color',C2(2,:))
hold on
plot(0,[per_ask_side_minus(10)],'^-.','color',C2(1,:))
lh = legend([h1 h2 h3 h4],'ASK immediate change (type A)','ASK immediate change (type B)','ASK reversion (type A)','ASK reversion (type B)');
legend boxoff
if lind == 1
set(lh, 'Position', [.5+.02 .885-0.277 .5 .01])
xlabel('\DeltaA_0')
ylabel('E[R_a(\tau)|\DeltaA_0],   E[\Deltaa_0|\DeltaA_0]')
text(-0.1,1.1,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlim([-5.5 5.5])
grid off
grid off
box on
end

subplot(2,2,3)
hold on
h1 = plot(imm_ask(2:end,1),imm_mid(2:end,1),'*-','color',C2(2,:),'LineWidth',lw);
hold on
plot(imm_ask(1:2,1),imm_mid(1:2,1),'-','color',C2(2,:),'LineWidth',lw);
hold on
plot(imm_ask(1,1),imm_mid(1,1),'d','color',C2(2,:),'LineWidth',lw);
hold on
plot(imm_ask(1:2,2),imm_mid(1:2,2),'--','color',C2(2,:),'LineWidth',lw)
hold on
h2 = plot(imm_ask(2:end,2),imm_mid(2:end,2),'o--','color',C2(2,:),'LineWidth',lw);
hold on
plot(imm_ask(2:end,1),per_mid(2:end,1),'*-','color',C2(1,:),'LineWidth',lw)
hold on
plot(imm_ask(1:2,1),per_mid(1:2,1),'-','color',C2(1,:),'LineWidth',lw)
hold on
plot(imm_ask(1,1),per_mid(1,1),'d','color',C2(1,:),'LineWidth',lw)
hold on
plot(imm_ask(1:2,2),per_mid(1:2,2),'--','color',C2(1,:),'LineWidth',lw)
hold on
plot(imm_ask(2:end,2),per_mid(2:end,2),'o--','color',C2(1,:),'LineWidth',lw)
hold on
plot(0:1,[imm_mid_side(10);imm_mid(1,1)],'-.','color',C2(2,:),'LineWidth',lw)
hold on
plot(0:1,[per_mid_side(10);per_mid(1,1)],'-.','color',C2(1,:),'LineWidth',lw)
hold on
plot(0,[imm_mid_side(10)],'^-.','color',C2(2,:))
hold on
plot(0,[per_mid_side(10)],'^-.','color',C2(1,:))
hold on

plot(imm_bid_minus(2:end,1),imm_mid_minus(2:end,1),'*-','color',C2(2,:),'LineWidth',lw)
hold on
plot(imm_bid_minus(1:2,1),imm_mid_minus(1:2,1),'-','color',C2(2,:),'LineWidth',lw)
hold on
plot(imm_bid_minus(1,1),imm_mid_minus(1,1),'d','color',C2(2,:),'LineWidth',lw)
hold on
plot(imm_bid_minus(1:2,2),imm_mid_minus(1:2,2),'--','color',C2(2,:),'LineWidth',lw)
hold on
plot(imm_bid_minus(2:end,2),imm_mid_minus(2:end,2),'o--','color',C2(2,:),'LineWidth',lw)
hold on
h3 = plot(imm_bid_minus(2:end,1),per_mid_minus(2:end,1),'*-','color',C2(1,:),'LineWidth',lw);
hold on
plot(imm_bid_minus(1:2,1),per_mid_minus(1:2,1),'-','color',C2(1,:),'LineWidth',lw);
hold on
plot(imm_bid_minus(1,1),per_mid_minus(1,1),'d','color',C2(1,:),'LineWidth',lw);
hold on
plot(imm_bid_minus(1:2,2),per_mid_minus(1:2,2),'--','color',C2(1,:),'LineWidth',lw)
hold on
h4 = plot(imm_bid_minus(2:end,2),per_mid_minus(2:end,2),'o--','color',C2(1,:),'LineWidth',lw);
hold on
plot([0 -1],[imm_mid_side_minus(10);imm_mid_minus(1,1)],'-.','color',C2(2,:),'LineWidth',lw)
hold on
plot([0 -1],[per_mid_side_minus(10);per_mid_minus(1,1)],'-.','color',C2(1,:),'LineWidth',lw)
hold on
plot(0,[imm_mid_side_minus(10)],'^-.','color',C2(2,:))
hold on
plot(0,[per_mid_side_minus(10)],'^-.','color',C2(1,:))
if lind == 1
lh = legend([h1 h2 h3 h4],'MID immediate change (type A)','MID immediate change (type B)','MID reversion (type A)','MID reversion (type B)');
legend boxoff
set(lh, 'Position', [.01 .385 .5 .01])
xlabel('\DeltaA_0')
ylabel('E[R_m(\tau)|\DeltaA_0],   E[\Deltam_0|\DeltaA_0]')
text(-0.1,1.1,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlim([-5.5 5.5])
grid off
grid off
box on
end
subplot(2,2,4)
hold on
h1 = plot(imm_ask(2:end,1),imm_spread(2:end,1),'*-','color',C2(2,:),'LineWidth',lw);
hold on
plot(imm_ask(1:2,1),imm_spread(1:2,1),'-','color',C2(2,:),'LineWidth',lw);
hold on
plot(imm_ask(1,1),imm_spread(1,1),'d','color',C2(2,:),'LineWidth',lw);
hold on
plot(imm_ask(1:2,2),imm_spread(1:2,2),'--','color',C2(2,:),'LineWidth',lw)
hold on
h2 = plot(imm_ask(2:end,2),imm_spread(2:end,2),'o--','color',C2(2,:),'LineWidth',lw);
hold on
plot(imm_ask(2:end,1),per_spread(2:end,1),'*-','color',C2(1,:),'LineWidth',lw)
hold on
plot(imm_ask(1:2,1),per_spread(1:2,1),'-','color',C2(1,:),'LineWidth',lw)
hold on
plot(imm_ask(1,1),per_spread(1,1),'d','color',C2(1,:),'LineWidth',lw)
hold on
plot(imm_ask(1:2,2),per_spread(1:2,2),'--','color',C2(1,:),'LineWidth',lw)
hold on
plot(imm_ask(2:end,2),per_spread(2:end,2),'o--','color',C2(1,:),'LineWidth',lw)
hold on
plot(0:1,[imm_spread_side(10);imm_spread(1,1)],'-.','color',C2(2,:),'LineWidth',lw)
hold on
plot(0:1,[per_spread_side(10);per_spread(1,1)],'-.','color',C2(1,:),'LineWidth',lw)
hold on
h5 = plot(0,[imm_spread_side(10)],'^-.','color',C2(2,:));
hold on
plot(0,[per_spread_side(10)],'^-.','color',C2(1,:))
hold on

plot(imm_bid_minus(2:end,1),imm_spread_minus(2:end,1),'*-','color',C2(2,:),'LineWidth',lw)
hold on
plot(imm_bid_minus(1:2,1),imm_spread_minus(1:2,1),'-','color',C2(2,:),'LineWidth',lw)
hold on
plot(imm_bid_minus(1,1),imm_spread_minus(1,1),'d','color',C2(2,:),'LineWidth',lw)
hold on
plot(imm_bid_minus(1:2,2),imm_spread_minus(1:2,2),'--','color',C2(2,:),'LineWidth',lw)
hold on
plot(imm_bid_minus(2:end,2),imm_spread_minus(2:end,2),'o--','color',C2(2,:),'LineWidth',lw)
hold on
h3 = plot(imm_bid_minus(2:end,1),per_spread_minus(2:end,1),'*-','color',C2(1,:),'LineWidth',lw);
hold on
plot(imm_bid_minus(1:2,1),per_spread_minus(1:2,1),'-','color',C2(1,:),'LineWidth',lw);
hold on
plot(imm_bid_minus(1,1),per_spread_minus(1,1),'d','color',C2(1,:),'LineWidth',lw);
hold on
plot(imm_bid_minus(1:2,2),per_spread_minus(1:2,2),'--','color',C2(1,:),'LineWidth',lw)
hold on
h4 = plot(imm_bid_minus(2:end,2),per_spread_minus(2:end,2),'o--','color',C2(1,:),'LineWidth',lw);
hold on
plot([0 -1],[imm_spread_side_minus(10);imm_spread_minus(1,1)],'-.','color',C2(2,:),'LineWidth',lw)
hold on
plot([0 -1],[per_spread_side_minus(10);per_spread_minus(1,1)],'-.','color',C2(1,:),'LineWidth',lw)
hold on
plot(0,[imm_spread_side_minus(10)],'^-.','color',C2(2,:))
hold on
plot(0,[per_spread_side_minus(10)],'^-.','color',C2(1,:))
lh = legend([h1 h2 h3 h4 h5],'SPREAD immediate change (type A)','SPREAD immediate change (type B)','SPREAD reversion (type A)','SPREAD reversion (type B)');
if lind == 1
legend boxoff
set(lh, 'Position', [.5 .385 .5 .01])
xlabel('\DeltaA_0')
ylabel('E[R_s(\tau)|\DeltaA_0],   E[\Deltas_0|\DeltaA_0]')
text(-0.1,1.1,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlim([-5.5 5.5])
grid off
grid off
box on


dx=0.06;
dy=0.06;
x=(1-4*dx)/1.9;
y=(1-4*dy)/2;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure2,'Type','axes');
set(AxesHandle(4),'Position',[dx,y+3*dy,x,y]);
set(AxesHandle(3),'Position',[x+2.5*dx,y+3*dy,x,y]);
set(AxesHandle(2),'Position',[dx,dy,x,y]);
set(AxesHandle(1),'Position',[x+2.5*dx,dy,x,y]);
set(gcf,'Position',[100, 100, 1000, 800])
end




% figure(3)
% % figure100 = figure('Position',[1 1 1200 500],'Color',[1 1 1]);
% set(gcf,'color','w')
% subplot(1,3,1)
% hold on
% plot(imm_ask(:,1),imm_ask(:,1)+per_ask(:,1),'b*-','LineWidth',lw)
% hold on
% plot(imm_ask(1:2,2),imm_ask(1:2,2)+per_ask(1:2,2),'b--','LineWidth',lw)
% hold on
% plot(imm_ask(2:end,2),imm_ask(2:end,2)+per_ask(2:end,2),'bo--','LineWidth',lw)
% hold on
% plot(imm_ask(:,1),imm_bid(:,1)+per_bid(:,1),'r*-','LineWidth',lw)
% hold on
% plot(imm_ask(1:2,2),imm_bid(1:2,2)+per_bid(1:2,2),'r--','LineWidth',lw)
% hold on
% plot(imm_ask(2:end,2),imm_bid(2:end,2)+per_bid(2:end,2),'ro--','LineWidth',lw)
% hold on
% plot(0:1,[imm_ask_side(10)+per_ask_side(10);imm_ask(1)+per_ask(1)],'b-.','LineWidth',lw)
% hold on
% plot(0:1,[imm_bid_side(10)+per_bid_side(10);imm_bid(1)+per_bid(1)],'r-.','LineWidth',lw)
% hold on
% plot(0,[imm_ask_side(10)+per_ask_side(10)],'b^-.')
% hold on
% plot(0,[imm_bid_side(10)+per_bid_side(10)],'r^-.')
% 
% xlim([0 5])
% legend('ASK permanent change (type A)','ASK permanent change (type B)','BID permanent change (type A)','BID permanent change (type B)')
% xlabel('\Deltaa_0')
% ylabel('E[P_a(\tau)|\Deltaa_0],  E[P_b(\tau)|\Deltaa_0]')
% text(-0.1,1.1,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
% 
% subplot(1,3,2)
% hold on
% plot(imm_ask(:,1),imm_mid(:,1)+per_mid(:,1),'k*-','LineWidth',lw)
% hold on
% plot(imm_ask(1:2,2),imm_mid(1:2,2)+per_mid(1:2,2),'k--','LineWidth',lw)
% hold on
% plot(imm_ask(2:end,2),imm_mid(2:end,2)+per_mid(2:end,2),'ko--','LineWidth',lw)
% hold on
% plot(0:1,[imm_mid_side(10)+per_mid_side(10);imm_mid(1)+per_mid(1)],'k-.','LineWidth',lw)
% hold on
% plot(0,[imm_mid_side(10)+per_mid_side(10)],'k^-.')
% xlim([0 5])
% legend('MID permanent change (type A)','MID permanent change (type B)')
% xlabel('\Deltaa_0')
% ylabel('E[P_M(\tau)|\Deltaa_0]')
% text(-0.1,1.1,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
% 
% subplot(1,3,3)
% hold on
% plot(imm_ask(:,1),imm_spread(:,1)+per_spread(:,1),'k*-','LineWidth',lw)
% hold on
% plot(imm_ask(1:2,2),imm_spread(1:2,2)+per_spread(1:2,2),'k--','LineWidth',lw)
% hold on
% plot(imm_ask(2:end,2),imm_spread(2:end,2)+per_spread(2:end,2),'ko--','LineWidth',lw)
% hold on
% plot(0:1,[imm_spread_side(10)+per_spread_side(10);imm_spread(1)+per_spread(1)],'k-.','LineWidth',lw)
% hold on
% plot(0,[imm_spread_side(10)+per_spread_side(10)],'k^-.')
% xlim([0 5])
% legend('SPREAD permanent change (type A)','SPREAD permanent change (type B)')
% xlabel('\Deltaa_0')
% ylabel('E[P_S(\tau)|\Deltaa_0]')
% text(-0.1,1.1,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top')
% 
% dx=0.05;
% dy=0.05;
% x=(1-5*dx)/2.9;
% y=(1-3*dy)/1.1;
% % AxesHandle=findobj(figure100,'Type','axes');
% % set(AxesHandle(1),'Position',[0.01+2*x+4*dx,0.08+dy,x,y]);
% % set(AxesHandle(2),'Position',[0.01+x+2.5*dx,0.08+dy,x,y]);
% % set(AxesHandle(3),'Position',[0.01+dx,0.08+dy,x,y]);


figure3 = figure(3);
% figure100 = figure('Position',[1 1 1200 500],'Color',[1 1 1]);
set(gcf,'color','w')
subplot(2,2,1)
hold on
plot(imm_ask(2:end,1),imm_bid(2:end,1)+per_bid(2:end,1),'*-','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_ask(2:end,2),imm_bid(2:end,2)+per_bid(2:end,2),'o--','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_ask(1:2,1),imm_bid(1:2,1)+per_bid(1:2,1),'-','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_ask(1,1),imm_bid(1,1)+per_bid(1,1),'d','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_bid_minus(2:end,1),abs(imm_bid_minus(2:end,1)+per_bid_minus(2:end,1)),'*-','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_bid_minus(1:2,1),abs(imm_bid_minus(1:2,1)+per_bid_minus(1:2,1)),'-','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_bid_minus(1,1),abs(imm_bid_minus(1,1)+per_bid_minus(1,1)),'d','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_bid_minus(2:end,2),abs(imm_bid_minus(2:end,2)+per_bid_minus(2:end,2)),'o--','LineWidth',lw,'color',[c1 c2 c3])

hold on
plot(imm_ask(1:2,2),imm_bid(1:2,2)+per_bid(1:2,2),'--','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(0:1,[imm_bid_side(10)+per_bid_side(10);imm_bid(1)+per_bid(1)],'-.','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(0,[imm_bid_side(10)+per_bid_side(10)],'^-.','color',[c1 c2 c3])

hold on
plot(imm_bid_minus(1:2,2),abs(imm_bid_minus(1:2,2)+per_bid_minus(1:2,2)),'--','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot([0,-1],[abs(imm_bid_side_minus(10)+per_bid_side_minus(10));abs(imm_bid_minus(1)+per_bid_minus(1))],'-.','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(0,[abs(imm_bid_side_minus(10)+per_bid_side_minus(10))],'^-.','color',[c1 c2 c3])
grid off
grid off
box on

if lind == 1
xlim([-5.5 5.5])
lh = legend('BID long-term change (type A)','BID long-term change (type B)');
legend boxoff
set(lh, 'Position', [.01 .9 .5 .01])
xlabel('\DeltaA_0')
ylabel('|E[I_b(\tau)|\DeltaA_0]|')
text(-0.1,1.1,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
end

subplot(2,2,2)
hold on
plot(imm_ask(2:end,1),imm_ask(2:end,1)+per_ask(2:end,1),'*-','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_ask(2:end,2),imm_ask(2:end,2)+per_ask(2:end,2),'o--','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_ask(1:2,1),imm_ask(1:2,1)+per_ask(1:2,1),'-','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_ask(1,1),imm_ask(1,1)+per_ask(1,1),'d','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_bid_minus(2:end,1),abs(imm_ask_minus(2:end,1)+per_ask_minus(2:end,1)),'*-','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_bid_minus(1:2,1),abs(imm_ask_minus(1:2,1)+per_ask_minus(1:2,1)),'-','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_bid_minus(1,1),abs(imm_ask_minus(1,1)+per_ask_minus(1,1)),'d','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_bid_minus(2:end,2),abs(imm_ask_minus(2:end,2)+per_ask_minus(2:end,2)),'o--','LineWidth',lw,'color',[c1 c2 c3])

hold on
plot(imm_ask(1:2,2),imm_ask(1:2,2)+per_ask(1:2,2),'--','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(0:1,[imm_ask_side(10)+per_ask_side(10);imm_ask(1)+per_ask(1)],'-.','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(0,[imm_ask_side(10)+per_ask_side(10)],'^-.','color',[c1 c2 c3])

hold on
plot(imm_bid_minus(1:2,2),abs(imm_ask_minus(1:2,2)+per_ask_minus(1:2,2)),'--','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot([0 -1],[abs(imm_ask_side_minus(10)+per_ask_side_minus(10));abs(imm_ask_minus(1)+per_ask_minus(1))],'-.','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(0,[abs(imm_ask_side_minus(10)+per_ask_side_minus(10))],'^-.','color',[c1 c2 c3])

if lind == 1
xlim([-5.5 5.5])
lh = legend('ASK long-term change (type A)','ASK long-term change (type B)');
legend boxoff
set(lh, 'Position', [.5 .9 .5 .01])
xlabel('\DeltaA_0')
ylabel('|E[I_a(\tau)|\DeltaA_0]|')
text(-0.1,1.1,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
grid off
grid off
box on
end

subplot(2,2,3)
hold on
plot(imm_ask(2:end,1),imm_mid(2:end,1)+per_mid(2:end,1),'*-','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_ask(2:end,2),imm_mid(2:end,2)+per_mid(2:end,2),'o--','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_ask(1:2,1),imm_mid(1:2,1)+per_mid(1:2,1),'-','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_ask(1,1),imm_mid(1,1)+per_mid(1,1),'d','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_bid_minus(2:end,1),abs(imm_mid_minus(2:end,1)+per_mid_minus(2:end,1)),'*-','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_bid_minus(1:2,1),abs(imm_mid_minus(1:2,1)+per_mid_minus(1:2,1)),'-','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_bid_minus(1,1),abs(imm_mid_minus(1,1)+per_mid_minus(1,1)),'d','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_bid_minus(2:end,2),abs(imm_mid_minus(2:end,2)+per_mid_minus(2:end,2)),'o--','LineWidth',lw,'color',[c1 c2 c3])

hold on
plot(imm_ask(1:2,2),imm_mid(1:2,2)+per_mid(1:2,2),'--','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(0:1,[imm_mid_side(10)+per_mid_side(10);imm_mid(1)+per_mid(1)],'-.','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(0,[imm_mid_side(10)+per_mid_side(10)],'^-.','color',[c1 c2 c3])

hold on
plot(imm_bid_minus(1:2,2),abs(imm_mid_minus(1:2,2)+per_mid_minus(1:2,2)),'--','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot([0 -1],[abs(imm_mid_side_minus(10)+per_mid_side_minus(10));abs(imm_mid_minus(1)+per_mid_minus(1))],'-.','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(0,abs(imm_mid_side_minus(10)+per_mid_side_minus(10)),'^-.','color',[c1 c2 c3])
xlim([-5.5 5.5])
lh = legend('MID long-term change (type A)','MID long-term change (type B)');
legend boxoff
if lind == 1
set(lh, 'Position', [.01 .4 .5 .01])
xlabel('\DeltaA_0')
ylabel('|E[I_M(\tau)|\DeltaA_0]|')
text(-0.1,1.1,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top')
grid off
grid off
box on
end

subplot(2,2,4)
hold on
plot(imm_ask(2:end,1),imm_spread(2:end,1)+per_spread(2:end,1),'*-','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_ask(2:end,2),imm_spread(2:end,2)+per_spread(2:end,2),'o--','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_ask(1:2,1),imm_spread(1:2,1)+per_spread(1:2,1),'-','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_ask(1,1),imm_spread(1,1)+per_spread(1,1),'d','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_bid_minus(2:end,1),imm_spread_minus(2:end,1)+per_spread_minus(2:end,1),'*-','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_bid_minus(1:2,1),imm_spread_minus(1:2,1)+per_spread_minus(1:2,1),'-','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_bid_minus(1,1),imm_spread_minus(1,1)+per_spread_minus(1,1),'d','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(imm_bid_minus(2:end,2),imm_spread_minus(2:end,2)+per_spread_minus(2:end,2),'o--','LineWidth',lw,'color',[c1 c2 c3])


hold on
plot(imm_ask(1:2,2),imm_spread(1:2,2)+per_spread(1:2,2),'--','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(0:1,[imm_spread_side(10)+per_spread_side(10);imm_spread(1)+per_spread(1)],'-.','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(0,[imm_spread_side(10)+per_spread_side(10)],'^-.','color',[c1 c2 c3])

hold on
plot(imm_bid_minus(1:2,2),imm_spread_minus(1:2,2)+per_spread_minus(1:2,2),'--','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot([0 -1],[imm_spread_side_minus(10)+per_spread_side_minus(10);imm_spread_minus(1)+per_spread_minus(1)],'-.','LineWidth',lw,'color',[c1 c2 c3])
hold on
plot(0,[imm_spread_side_minus(10)+per_spread_side_minus(10)],'^-.','color',[c1 c2 c3])
xlim([-5.5 5.5])
lh = legend('SPREAD long-term change (type A)','SPREAD long-term change (type B)');
legend boxoff
if lind == 1
set(lh, 'Position', [.5 .4 .5 .01])
xlabel('\DeltaA_0')
ylabel('E[I_S(\tau)|\DeltaA_0]')
text(-0.1,1.1,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top')
grid off
grid off
box on

dx=0.06;
dy=0.06;
x=(1-4*dx)/1.9;
y=(1-4*dy)/2;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure3,'Type','axes');
set(AxesHandle(4),'Position',[dx,y+3*dy,x,y]);
set(AxesHandle(3),'Position',[x+2.5*dx,y+3*dy,x,y]);
set(AxesHandle(2),'Position',[dx,dy,x,y]);
set(AxesHandle(1),'Position',[x+2.5*dx,dy,x,y]);
set(gcf,'Position',[100, 100, 1000, 800])
end





figure4 = figure(4);
set(gcf,'color','w')
if day_start == 1
    hold on
    f1 = plot(imm_ask(2:end,1),imm_spread(2:end,1)+per_spread(2:end,1),'*-','LineWidth',lw,'color',C4(ci,:));
    hold on
    plot(imm_ask(1:2,1),imm_spread(1:2,1)+per_spread(1:2,1),'-','LineWidth',lw,'color',C4(ci,:));
    hold on
    plot(imm_ask(1,1),imm_spread(1,1)+per_spread(1,1),'d','LineWidth',lw,'color',C4(ci,:));
    hold on
    f2 = plot(imm_ask(2:end,2),imm_spread(2:end,2)+per_spread(2:end,2),'o--','LineWidth',lw,'color',C4(ci,:));
    hold on
    plot(imm_bid_minus(2:end,1),imm_spread_minus(2:end,1)+per_spread_minus(2:end,1),'*-','LineWidth',lw,'color',C4(ci,:));
    hold on
    plot(imm_bid_minus(1:2,1),imm_spread_minus(1:2,1)+per_spread_minus(1:2,1),'-','LineWidth',lw,'color',C4(ci,:));
    hold on
    plot(imm_bid_minus(1,1),imm_spread_minus(1,1)+per_spread_minus(1,1),'d','LineWidth',lw,'color',C4(ci,:));
    hold on
    plot(imm_bid_minus(2:end,2),imm_spread_minus(2:end,2)+per_spread_minus(2:end,2),'o--','LineWidth',lw,'color',C4(ci,:));

    hold on
    plot(imm_ask(1:2,2),imm_spread(1:2,2)+per_spread(1:2,2),'--','LineWidth',lw,'color',C4(ci,:))
    hold on
    plot(0:1,[imm_spread_side(10)+per_spread_side(10);imm_spread(1)+per_spread(1)],'-.','LineWidth',lw,'color',C4(ci,:))
    hold on
    plot(0,[imm_spread_side(10)+per_spread_side(10)],'^-.','color',C4(ci,:))

    hold on
    plot(imm_bid_minus(1:2,2),imm_spread_minus(1:2,2)+per_spread_minus(1:2,2),'--','LineWidth',lw,'color',C4(ci,:))
    hold on
    plot([0 -1],[imm_spread_side_minus(10)+per_spread_side_minus(10);imm_spread_minus(1)+per_spread_minus(1)],'-.','LineWidth',lw,'color',C4(ci,:))
    hold on
    plot(0,[imm_spread_side_minus(10)+per_spread_side_minus(10)],'^-.','color',C4(ci,:))
    hold on
elseif day_start == 32
    hold on
    f3 = plot(imm_ask(2:end,1),imm_spread(2:end,1)+per_spread(2:end,1),'*-','LineWidth',lw,'color',C4(ci,:));
    hold on
    plot(imm_ask(1:2,1),imm_spread(1:2,1)+per_spread(1:2,1),'-','LineWidth',lw,'color',C4(ci,:));
    hold on
    plot(imm_ask(1,1),imm_spread(1,1)+per_spread(1,1),'d','LineWidth',lw,'color',C4(ci,:));
    hold on
    f4 = plot(imm_ask(2:end,2),imm_spread(2:end,2)+per_spread(2:end,2),'o--','LineWidth',lw,'color',C4(ci,:));
    hold on
    plot(imm_bid_minus(2:end,1),imm_spread_minus(2:end,1)+per_spread_minus(2:end,1),'*-','LineWidth',lw,'color',C4(ci,:))
    hold on
    plot(imm_bid_minus(1:2,1),imm_spread_minus(1:2,1)+per_spread_minus(1:2,1),'-','LineWidth',lw,'color',C4(ci,:))
    hold on
    plot(imm_bid_minus(1,1),imm_spread_minus(1,1)+per_spread_minus(1,1),'d','LineWidth',lw,'color',C4(ci,:))
    hold on
    plot(imm_bid_minus(2:end,2),imm_spread_minus(2:end,2)+per_spread_minus(2:end,2),'o--','LineWidth',lw,'color',C4(ci,:))

    hold on
    plot(imm_ask(1:2,2),imm_spread(1:2,2)+per_spread(1:2,2),'--','LineWidth',lw,'color',C4(ci,:))
    hold on
    plot(0:1,[imm_spread_side(10)+per_spread_side(10);imm_spread(1)+per_spread(1)],'-.','LineWidth',lw,'color',C4(ci,:))
    hold on
    plot(0,[imm_spread_side(10)+per_spread_side(10)],'^-.','color',C4(ci,:))

    hold on
    plot(imm_bid_minus(1:2,2),imm_spread_minus(1:2,2)+per_spread_minus(1:2,2),'--','LineWidth',lw,'color',C4(ci,:))
    hold on
    plot([0 -1],[imm_spread_side_minus(10)+per_spread_side_minus(10);imm_spread_minus(1)+per_spread_minus(1)],'-.','LineWidth',lw,'color',C4(ci,:))
    hold on
    plot(0,[imm_spread_side_minus(10)+per_spread_side_minus(10)],'^-.','color',C4(ci,:))
    hold on
elseif day_start == 140
    hold on
    f5 = plot(imm_ask(2:end,1),imm_spread(2:end,1)+per_spread(2:end,1),'*-','LineWidth',lw,'color',C4(ci,:));
    hold on
    plot(imm_ask(1:2,1),imm_spread(1:2,1)+per_spread(1:2,1),'-','LineWidth',lw,'color',C4(ci,:));
    hold on
    plot(imm_ask(1,1),imm_spread(1,1)+per_spread(1,1),'d','LineWidth',lw,'color',C4(ci,:));
    hold on
    f6 = plot(imm_ask(2:end,2),imm_spread(2:end,2)+per_spread(2:end,2),'o--','LineWidth',lw,'color',C4(ci,:));
    hold on
    plot(imm_bid_minus(2:end,1),imm_spread_minus(2:end,1)+per_spread_minus(2:end,1),'*-','LineWidth',lw,'color',C4(ci,:))
    hold on
    plot(imm_bid_minus(1:2,1),imm_spread_minus(1:2,1)+per_spread_minus(1:2,1),'-','LineWidth',lw,'color',C4(ci,:))
    hold on
    plot(imm_bid_minus(1,1),imm_spread_minus(1,1)+per_spread_minus(1,1),'d','LineWidth',lw,'color',C4(ci,:))
    hold on
    plot(imm_bid_minus(2:end,2),imm_spread_minus(2:end,2)+per_spread_minus(2:end,2),'o--','LineWidth',lw,'color',C4(ci,:))

    hold on
    plot(imm_ask(1:2,2),imm_spread(1:2,2)+per_spread(1:2,2),'--','LineWidth',lw,'color',C4(ci,:))
    hold on
    plot(0:1,[imm_spread_side(10)+per_spread_side(10);imm_spread(1)+per_spread(1)],'-.','LineWidth',lw,'color',C4(ci,:))
    hold on
    plot(0,[imm_spread_side(10)+per_spread_side(10)],'^-.','color',C4(ci,:))

    hold on
    plot(imm_bid_minus(1:2,2),imm_spread_minus(1:2,2)+per_spread_minus(1:2,2),'--','LineWidth',lw,'color',C4(ci,:))
    hold on
    plot([0 -1],[imm_spread_side_minus(10)+per_spread_side_minus(10);imm_spread_minus(1)+per_spread_minus(1)],'-.','LineWidth',lw,'color',C4(ci,:))
    hold on
    plot(0,[imm_spread_side_minus(10)+per_spread_side_minus(10)],'^-.','color',C4(ci,:))
    hold on
end

if lind == 1
xlim([-5.5 5.5])
legend([f1,f2,f3,f4,f5,f6],{'SPREAD (Type A, BC)','SPREAD (Type B, BC)','SPREAD (Type A, DC)','SPREAD (Type B, DC)','SPREAD (Type A, AC)','SPREAD (Type B, AC)'})
%lh = legend('SPREAD permanent change (type A)','SPREAD permanent change (type B)');
legend boxoff
set(lh, 'Position', [.5 .4 .5 .01])
xlabel('\DeltaA_0')
ylabel('E[I_S(\tau)|\DeltaA_0]')
% text(-0.1,1.1,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top')
grid off
grid off
box on

dx=0.06;
dy=0.06;
x=(1-4*dx)/0.9;
y=(1-4*dy)/0.0;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure4,'Type','axes');
% set(AxesHandle(1),'Position',[dx,dy,x,y]);
set(gcf,'Position',[100, 100, 500, 370])
end


% 
% figure
% set(gcf,'color','w')
% subplot(1,3,1)
% plot(imm_ask(:,1),imm_ask_vol_result(:,1),'b*-')
% hold on
% plot(imm_ask(:,2),imm_ask_vol_result(:,2),'bo--')
% hold on
% plot(imm_ask(:,1),imm_bid_vol_result(:,1),'r*-')
% hold on
% plot(imm_ask(:,2),imm_bid_vol_result(:,2),'ro--')
% legend('ASK VOL immediate change (type A)','ASK VOL immediate change (type B)','BID VOL immediate change (type A)','BID VOL immediate change (type B)')
% xlabel('\Deltaa_0')
% ylabel('E[R(\tau)|\Deltaa_0],   E[IM(0)|\Deltaa_0]','FontSize',16)
% 
% subplot(1,3,2)
% plot(imm_ask(:,1),imm_ask_gap_result(:,1),'b*-')
% hold on
% plot(imm_ask(:,2),imm_ask_gap_result(:,2),'bo--')
% hold on
% plot(imm_ask(:,1),imm_bid_gap_result(:,1),'r*-')
% hold on
% plot(imm_ask(:,2),imm_bid_gap_result(:,2),'ro--')
% legend('ASK GAP immediate change (type A)','ASK GAP immediate change (type B)','BID GAP immediate change (type A)','BID GAP immediate change (type B)')
% xlabel('\Deltaa_0')
% ylabel('E[R(\tau)|\Deltaa_0],   E[IM(0)|\Deltaa_0]','FontSize',16)
% 
% subplot(1,3,3)
% plot(imm_ask(:,1),imm_ask_shock_result(:,1),'b*-')
% hold on
% plot(imm_ask(:,2),imm_ask_shock_result(:,2),'bo--')
% hold on
% plot(imm_ask(:,1),imm_bid_shock_result(:,1),'r*-')
% hold on
% plot(imm_ask(:,2),imm_bid_shock_result(:,2),'ro--')
% legend('ASK SHOCK immediate change (type A)','ASK SHOCK immediate change (type B)','BID SHOCK immediate change (type A)','BID SHOCK immediate change (type B)')
% xlabel('\Deltaa_0')
% ylabel('E[R(\tau)|\Deltaa_0],   E[IM(0)|\Deltaa_0]','FontSize',16)
% 
% 
% 
% figure
% set(gcf,'color','w')
% subplot(1,3,1)
% plot(imm_ask(:,1),per_ask_vol_result(:,1),'b*-')
% hold on
% plot(imm_ask(:,2),per_ask_vol_result(:,2),'bo--')
% hold on
% plot(imm_ask(:,1),per_bid_vol_result(:,1),'r*-')
% hold on
% plot(imm_ask(:,2),per_bid_vol_result(:,2),'ro--')
% hold on
% legend('ASK VOL reversion change (type A)','ASK VOL reversion change (type B)','BID VOL reversion change (type A)','BID VOL reversion change (type B)')
% xlabel('\Deltaa_0')
% ylabel('E[R(\tau)|\Deltaa_0],   E[IM(0)|\Deltaa_0]','FontSize',16)
% 
% subplot(1,3,2)
% plot(imm_ask(:,1),per_ask_gap_result(:,1),'b*-')
% hold on
% plot(imm_ask(:,2),per_ask_gap_result(:,2),'bo--')
% hold on
% plot(imm_ask(:,1),per_bid_gap_result(:,1),'r*-')
% hold on
% plot(imm_ask(:,2),per_bid_gap_result(:,2),'ro--')
% legend('ASK GAP reversion change (type A)','ASK GAP reversion change (type B)','BID GAP reversion change (type A)','BID GAP reversion change (type B)')
% xlabel('\Deltaa_0')
% ylabel('E[R(\tau)|\Deltaa_0],   E[IM(0)|\Deltaa_0]','FontSize',16)
% 
% subplot(1,3,3)
% plot(imm_ask(:,1),per_ask_shock_result(:,1),'b*-')
% hold on
% plot(imm_ask(:,2),per_ask_shock_result(:,2),'bo--')
% hold on
% plot(imm_ask(:,1),per_bid_shock_result(:,1),'r*-')
% hold on
% plot(imm_ask(:,2),per_bid_shock_result(:,2),'ro--')
% legend('ASK SHOCK reversion change (type A)','ASK SHOCK reversion change (type B)','BID SHOCK reversion change (type A)','BID SHOCK reversion change (type B)')
% xlabel('\Deltaa_0')
% ylabel('E[R(\tau)|\Deltaa_0],   E[IM(0)|\Deltaa_0]','FontSize',16)
% 
% 

% figure
% set(gcf,'color','w')
% subplot(1,3,1)
% plot(imm_ask(:,1),imm_ask_vol_result(:,1)+per_ask_vol_result(:,1),'b*-')
% hold on
% plot(imm_ask(:,2),imm_ask_vol_result(:,2)+per_ask_vol_result(:,2),'bo--')
% hold on
% plot(imm_ask(:,1),imm_bid_vol_result(:,1)+per_bid_vol_result(:,1),'r*-')
% hold on
% plot(imm_ask(:,2),imm_bid_vol_result(:,2)+per_bid_vol_result(:,2),'ro--')
% hold on
% legend('ASK VOL permanent change (type A)','ASK VOL permanent change (type B)','BID VOL permanent change (type A)','BID VOL permanent change (type B)')
% xlabel('\Deltaa_0')
% ylabel('E[R(\tau)|\Deltaa_0],   E[IM(0)|\Deltaa_0]','FontSize',16)
% 
% subplot(1,3,2)
% plot(imm_ask(:,1),imm_ask_gap_result(:,1)+per_ask_gap_result(:,1),'b*-')
% hold on
% plot(imm_ask(:,2),imm_ask_gap_result(:,2)+per_ask_gap_result(:,2),'bo--')
% hold on
% plot(imm_ask(:,1),imm_bid_gap_result(:,1)+per_bid_gap_result(:,1),'r*-')
% hold on
% plot(imm_ask(:,2),imm_bid_gap_result(:,2)+per_bid_gap_result(:,2),'ro--')
% legend('ASK GAP permanent change (type A)','ASK GAP permanent change (type B)','BID GAP permanent change (type A)','BID GAP permanent change (type B)')
% xlabel('\Deltaa_0')
% ylabel('E[R(\tau)|\Deltaa_0],   E[IM(0)|\Deltaa_0]','FontSize',16)
% 
% subplot(1,3,3)
% plot(imm_ask(:,1),imm_ask_shock_result(:,1)+per_ask_shock_result(:,1),'b*-')
% hold on
% plot(imm_ask(:,2),imm_ask_shock_result(:,2)+per_ask_shock_result(:,2),'bo--')
% hold on
% plot(imm_ask(:,1),imm_bid_shock_result(:,1)+per_bid_shock_result(:,1),'r*-')
% hold on
% plot(imm_ask(:,2),imm_bid_shock_result(:,2)+per_bid_shock_result(:,2),'ro--')
% legend('ASK SHOCK permanent change (type A)','ASK SHOCK permanent change (type B)','BID SHOCK permanent change (type A)','BID SHOCK permanent change (type B)')
% xlabel('\Deltaa_0')
% ylabel('E[R(\tau)|\Deltaa_0],   E[IM(0)|\Deltaa_0]','FontSize',16)
% % % %% bid, ask price (1 figure) mid price
% % % iter=iter_set;
% % % f1=figure;
% % % set(gcf,'color','w')
% % % % set(f1,'Visible','off')
% % % % set(f1, 'PaperUnits', 'inches');
% % % % x_width=14 ;y_width=7;
% % % % set(f1, 'PaperPosition', [0 0 x_width y_width]);
% % % 
% % % 
% % % for i=1:NI
% % %     if i==1
% % %         ind = 1;
% % %     elseif i>1 && i<4
% % %         ind = 2;
% % %         col = (i-1)/3;
% % %     elseif i>3 && i<7
% % %         ind = 3;
% % %         col = (i-3)/4;
% % %     elseif i>6 && i<11
% % %         ind = 4;
% % %         col = (i-6)/5;
% % %     elseif i>10
% % %         ind = 5;
% % %         col = (i-10)/6;
% % %     end
% % %     
% % %     if i==1 || i==2 || i==4 || i==7 || i==11
% % %         flag_g = 1;
% % %     else
% % %         flag_g = 2;
% % %     end
% % %     
% % %     refer_use =  (bid_price_plus_ta{iter,1}{i,1}(t+1) + ask_price_plus_ta{iter,1}{i,1}(t+1))/2;
% % % %     refer_use = 0;
% % %     
% % %     subplot(1,5,ind)
% % % 
% % %     if flag_g==1
% % %         plot(-t:t,(bid_price_plus_ta{iter,1}{i,1} + ask_price_plus_ta{iter,1}{i,1})/2 - refer_use,'k-','LineWidth',1)
% % %         hold on
% % %         plot(0,(bid_price_plus_ta{iter,1}{i,1}(t+1) + ask_price_plus_ta{iter,1}{i,1}(t+1))/2 - refer_use,'k*','MarkerSize',10)
% % %         hold on
% % %         plot(0,(bid_price_plus_ta{iter,1}{i,1}(t+2) + ask_price_plus_ta{iter,1}{i,1}(t+2))/2 - refer_use,'k*','MarkerSize',10)
% % %         hold on
% % %     elseif flag_g==2
% % %         plot(-t:t,(bid_price_plus_ta{iter,1}{i,1} + ask_price_plus_ta{iter,1}{i,1})/2 - refer_use,'Color',[1 col col],'LineWidth',1)
% % %         hold on
% % %         plot(0,(bid_price_plus_ta{iter,1}{i,1}(t+1) + ask_price_plus_ta{iter,1}{i,1}(t+1))/2 - refer_use,'o','Color',[1 col col],'MarkerSize',10)
% % %         hold on
% % %         plot(0,(bid_price_plus_ta{iter,1}{i,1}(t+2) + ask_price_plus_ta{iter,1}{i,1}(t+2))/2 - refer_use,'o','Color',[1 col col],'MarkerSize',10)
% % %         hold on
% % %     end
% % %     ylim([-.5 5])
% % % %         set(gca,'xscale','log')
% % % %         set(gca,'yscale','log')
% % % %         hold on
% % % %     hold on
% % % %     plot([-t:1:-1,1:t],spread_minus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_minus_num{iter,1}(i,1),'b.')
% % % %     legend('minus initiated market order','plus initiated market order')
% % %     xlabel('lag')
% % %     ylabel('Price')
% % % end
% % % 
% % % 
% % % 
% % % % iter=iter_set;
% % % % f2=figure;
% % % % set(gcf,'color','w')
% % % % % set(f1,'Visible','off')
% % % % % set(f1, 'PaperUnits', 'inches');
% % % % % x_width=14 ;y_width=7;
% % % % % set(f1, 'PaperPosition', [0 0 x_width y_width]);
% % % % 
% % % % 
% % % % for i=1:NI
% % % %     if i>5 && i<10
% % % %         ind = i-4;
% % % %         flag = 1;
% % % %         flag_g = 1;
% % % %     elseif i>9
% % % %         ind = i-8;
% % % %         flag = 1;
% % % %         flag_g = 2;
% % % %     elseif i==1
% % % %         ind = 1;
% % % %         flag = 1;
% % % %         flag_g = 1;
% % % %     elseif i>1 && i<6
% % % %         flag = 0;
% % % %     end
% % % %     
% % % %     
% % % %     refer_use =  (bid_price_minus_ta{iter,1}{i,1}(t+1) + ask_price_minus_ta{iter,1}{i,1}(t+1))/2;
% % % % %     refer_use = 0;
% % % %     
% % % %     if flag==1
% % % %         subplot(1,5,ind)
% % % %         
% % % %         if flag_g==1
% % % %             plot(-t:t,(bid_price_minus_ta{iter,1}{i,1} + ask_price_minus_ta{iter,1}{i,1})/2 - refer_use,'k-','LineWidth',1)
% % % %             hold on
% % % %             plot(0,(bid_price_minus_ta{iter,1}{i,1}(t+1) + ask_price_minus_ta{iter,1}{i,1}(t+1))/2 - refer_use,'k*','MarkerSize',10)
% % % %             hold on
% % % %         elseif flag_g==2
% % % %             plot(-t:t,(bid_price_minus_ta{iter,1}{i,1} + ask_price_minus_ta{iter,1}{i,1})/2 - refer_use,'k.','LineWidth',1)
% % % %             hold on
% % % %             plot(0,(bid_price_minus_ta{iter,1}{i,1}(t+1) + ask_price_minus_ta{iter,1}{i,1}(t+1))/2 - refer_use,'ko','MarkerSize',10)
% % % %             hold on
% % % %         end
% % % %         ylim([-5.5 .5])
% % % % %         set(gca,'xscale','log')
% % % % %         set(gca,'yscale','log')
% % % % %         hold on
% % % %     %     hold on
% % % %     %     plot([-t:1:-1,1:t],spread_minus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_minus_num{iter,1}(i,1),'b.')
% % % %     %     legend('minus initiated market order','plus initiated market order')
% % % %         xlabel('lag')
% % % %         ylabel('Price')
% % % %     end
% % % % end
% % % 
% % % 
% % % 
% % % %% flow
% % % iter=iter_set;
% % % for ii=5:5
% % %     ind=ii+6;
% % %     f1=figure;
% % %     set(gcf,'color','w')
% % % %     set(f1,'Visible','off')
% % % %     set(f1, 'PaperUnits', 'inches');
% % % %     x_width=14 ;y_width=7;
% % % %     set(f1, 'PaperPosition', [0 0 x_width y_width]);
% % %     for i=1:15
% % %         subplot(3,5,i)
% % %         plot([-t:1:-1,1:t],flow_plus_ta{iter,1}{i,ind}([1:t,t+2:2*t+1],1),'r.','LineWidth',2) % flow when minus initiated market orders come
% % %         hold on
% % %         plot([-t:1:-1,1:t],flow_minus_ta{iter,1}{i,ind}([1:t,t+2:2*t+1],1),'b.','LineWidth',2) % flow when plus initiated market orders come
% % %         hold on
% % %         plot([0 0],[min(min(flow_plus_ta{iter,1}{i,ind}([1:t,t+2:2*t+1],1)),min(flow_minus_ta{iter,1}{i,ind}([1:t,t+2:2*t+1],1))) max(max(flow_plus_ta{iter,1}{i,ind}([1:t,t+2:2*t+1],1)),max(flow_minus_ta{iter,1}{i,ind}([1:t,t+2:2*t+1],1)))],'k','LineWidth',.5)
% % %         xlim([-50 50])
% % %         
% % % %         plot(0,flow_plus_total{iter,1}{i,ind}(t+1)./flow_plus_num{iter,1}(i,ind),'r*','MarkerSize',10)
% % % %         hold on
% % % %         plot(0,flow_minus_total{iter,1}{i,ind}(t+1)./flow_minus_num{iter,1}(i,ind),'b*','MarkerSize',10)
% % % %         hold on
% % % %         plot([-t:-1,1:t],flow_plus_total{iter,1}{i,ind}([1:t,t+2:2*t+1])./flow_plus_num{iter,1}(i,ind),'ro')
% % % %         hold on
% % % %         plot([-t:-1,1:t],flow_minus_total{iter,1}{i,ind}([1:t,t+2:2*t+1])./flow_minus_num{iter,1}(i,ind),'bo')
% % %     %     legend('minus initiated market order','plus initiated market order')
% % %         xlabel('lag')
% % %         ylabel('flow')
% % %     %     set(gca,'xscale','log')
% % %     %     set(gca,'yscale','log')
% % %     %     xlim([0 50])
% % % %         if i==1
% % % %             title('Type A, Market order')
% % % %         elseif i==2
% % % %             title('Type B, Market order')
% % % %         elseif i==3
% % % %             title('Type C, Market order')
% % % %         elseif i==4
% % % %             title('Type D, Market order')
% % % %         elseif i==5
% % % %             title('Type E, Market order')
% % % %         elseif i==6
% % % %             title('Type F, Market order')
% % % %         elseif i==7
% % % %             title('Type G, Market order')
% % % %         end
% % %     end
% % % %     saveas(f1,sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/f1/%d_%d_%d_%d_flow.png',firm,iter,id,ii))
% % % end
% % % 
% % % 
% % % 
% % % %% order book depth
% % % iter=iter_set;
% % % % ii=1;
% % % for ii=1:1
% % % %     f1=figure;
% % %     figure
% % %     set(gcf,'color','w')
% % % %     set(f1,'Visible','off')
% % % %     set(f1, 'PaperUnits', 'inches');
% % % %     x_width=14 ;y_width=7;
% % % %     set(f1, 'PaperPosition', [0 0 x_width y_width]);
% % %     for i=1:NI
% % %         subplot(3,5,i)
% % %         plot(-t:t,bid_vol_plus_ta{iter,1}{i,ii},'r.','LineWidth',2)
% % %         hold on
% % %         plot(-t:t,ask_vol_plus_ta{iter,1}{i,ii},'b.','LineWidth',2)
% % %         hold on
% % %         plot(0,bid_vol_plus_ta{iter,1}{i,ii}(t+1),'r*','MarkerSize',10)
% % %         hold on
% % %         plot(0,ask_vol_plus_ta{iter,1}{i,ii}(t+1),'b*','MarkerSize',10)
% % % %         hold on
% % % %         plot([-t:-1,1:t],bid_vol_plus_total{iter,1}{i,ii}([1:t,t+2:2*t+1])./bid_vol_plus_num{iter,1}(i,ii),'ro')
% % % %         hold on
% % % %         plot([-t:-1,1:t],ask_vol_plus_total{iter,1}{i,ii}([1:t,t+2:2*t+1])./ask_vol_plus_num{iter,1}(i,ii),'bo')
% % %     %     legend('best bid depth','best ask depth')
% % %         xlabel('lag')
% % %         ylabel('Volume')
% % % %         if i==1
% % % %             title('Type A, minus initiated market order')
% % % %         elseif i==2
% % % %             title('Type B, minus initiated market order')
% % % %         elseif i==3
% % % %             title('Type C, minus initiated market order')
% % % %         elseif i==4
% % % %             title('Type D, minus initiated market order')
% % % %         elseif i==5
% % % %             title('Type E, minus initiated market order')
% % % %         elseif i==6
% % % %             title('Type F, minus initiated market order')
% % % %         elseif i==7
% % % %             title('Type G, minus initiated market order')
% % % %         end
% % %     end
% % % %     saveas(f1,sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/f1/%d_%d_%d_%d_depth_plusside.png',firm,iter,id,ii))
% % % 
% % % % % %     f1=figure;
% % % % % %     set(gcf,'color','w')
% % % % % % %     set(f1,'Visible','off')
% % % % % % %     set(f1, 'PaperUnits', 'inches');
% % % % % % %     x_width=14 ;y_width=7;
% % % % % % %     set(f1, 'PaperPosition', [0 0 x_width y_width]);
% % % % % %     for i=1:NI-1
% % % % % %         subplot(4,4,i)
% % % % % %         plot(-t:t,bid_vol_minus_total{iter,1}{i,ii}./bid_vol_minus_num{iter,1}(i,ii),'r','LineWidth',2)
% % % % % %         hold on
% % % % % %         plot(-t:t,ask_vol_minus_total{iter,1}{i,ii}./ask_vol_minus_num{iter,1}(i,ii),'b','LineWidth',2)
% % % % % %         hold on
% % % % % %         plot(0,bid_vol_minus_total{iter,1}{i,ii}(t+1)./bid_vol_minus_num{iter,1}(i,ii),'r*','MarkerSize',10)
% % % % % %         hold on
% % % % % %         plot(0,ask_vol_minus_total{iter,1}{i,ii}(t+1)./ask_vol_minus_num{iter,1}(i,ii),'b*','MarkerSize',10)
% % % % % %         hold on
% % % % % %         plot([-t:-1,1:t],bid_vol_minus_total{iter,1}{i,ii}([1:t,t+2:2*t+1])./bid_vol_minus_num{iter,1}(i,ii),'ro')
% % % % % %         hold on
% % % % % %         plot([-t:-1,1:t],ask_vol_minus_total{iter,1}{i,ii}([1:t,t+2:2*t+1])./ask_vol_minus_num{iter,1}(i,ii),'bo')
% % % % % %     %     legend('best bid depth','best ask depth')
% % % % % %         xlabel('lag')
% % % % % %         ylabel('Volume')
% % % % % %         if i==1
% % % % % %             title('Type A, plus initiated market order')
% % % % % %         elseif i==2
% % % % % %             title('Type B, plus initiated market order')
% % % % % %         elseif i==3
% % % % % %             title('Type C, plus initiated market order')
% % % % % %         elseif i==4
% % % % % %             title('Type D, plus initiated market order')
% % % % % %         elseif i==5
% % % % % %             title('Type E, plus initiated market order')
% % % % % %         elseif i==6
% % % % % %             title('Type F, plus initiated market order')
% % % % % %         elseif i==7
% % % % % %             title('Type G, plus initiated market order')
% % % % % %         end
% % % % % %     end
% % % %     saveas(f1,sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/f1/%d_%d_%d_%d_depth_minusside.png',firm,iter,id,ii))
% % % end
% % % 
% % % 
% % % %% order book depth (layer ratio)
% % % iter=iter_set;
% % % % ii=1;
% % % for ii=1:1
% % % %     f1=figure;
% % %     figure
% % %     set(gcf,'color','w')
% % % %     set(f1,'Visible','off')
% % % %     set(f1, 'PaperUnits', 'inches');
% % % %     x_width=14 ;y_width=7;
% % % %     set(f1, 'PaperPosition', [0 0 x_width y_width]);
% % %     for i=1:NI-1
% % %         subplot(4,4,i)
% % %         plot(-t:t,(bid_vol_plus_total{iter,1}{i,ii}/bid_vol_plus_num{iter,1}(i,ii))./(bid_vol_plus_total{iter,1}{i,ii+1}/bid_vol_plus_num{iter,1}(i,ii+1)),'r.','LineWidth',2)
% % %         hold on
% % %         plot(-t:t,(ask_vol_plus_total{iter,1}{i,ii}/ask_vol_plus_num{iter,1}(i,ii))./(ask_vol_plus_total{iter,1}{i,ii+1}/ask_vol_plus_num{iter,1}(i,ii+1)),'b.','LineWidth',2)
% % %         hold on
% % % %         plot(0,bid_vol_plus_total{iter,1}{i,ii}(t+1)/bid_vol_plus_num{iter,1}(i,ii),'r*','MarkerSize',10)
% % % %         hold on
% % % %         plot(0,ask_vol_plus_total{iter,1}{i,ii}(t+1)/ask_vol_plus_num{iter,1}(i,ii),'b*','MarkerSize',10)
% % % %         hold on
% % % %         plot([-t:-1,1:t],bid_vol_plus_total{iter,1}{i,ii}([1:t,t+2:2*t+1])./bid_vol_plus_num{iter,1}(i,ii),'ro')
% % % %         hold on
% % % %         plot([-t:-1,1:t],ask_vol_plus_total{iter,1}{i,ii}([1:t,t+2:2*t+1])./ask_vol_plus_num{iter,1}(i,ii),'bo')
% % %     %     legend('best bid depth','best ask depth')
% % %         xlabel('lag')
% % %         ylabel('Ratio')
% % % %         if i==1
% % % %             title('Type A, minus initiated market order')
% % % %         elseif i==2
% % % %             title('Type B, minus initiated market order')
% % % %         elseif i==3
% % % %             title('Type C, minus initiated market order')
% % % %         elseif i==4
% % % %             title('Type D, minus initiated market order')
% % % %         elseif i==5
% % % %             title('Type E, minus initiated market order')
% % % %         elseif i==6
% % % %             title('Type F, minus initiated market order')
% % % %         elseif i==7
% % % %             title('Type G, minus initiated market order')
% % % %         end
% % %     end
% % % end
% % % 
% % % 
% % % %% order book depth
% % % iter=iter_set;
% % % % ii=1;
% % % for ii=1:1
% % %     f1=figure;
% % %     set(gcf,'color','w')
% % % %     set(f1,'Visible','off')
% % % %     set(f1, 'PaperUnits', 'inches');
% % % %     x_width=14 ;y_width=7;
% % % %     set(f1, 'PaperPosition', [0 0 x_width y_width]);
% % %     for iii=1:3
% % %         if iii==1
% % %             i=1;
% % %         elseif iii==2
% % %             i=7;
% % %         elseif iii==3
% % %             i=11;
% % %         end
% % %         subplot(1,3,iii)
% % %         plot(-t:t,bid_vol_plus_total{iter,1}{i,ii}./bid_vol_plus_num{iter,1}(i,ii),'k','LineWidth',2)
% % %         hold on
% % %         plot(-t:t,ask_vol_plus_total{iter,1}{i,ii}./ask_vol_plus_num{iter,1}(i,ii),'k--','LineWidth',2)
% % %         hold on
% % %         plot(0,bid_vol_plus_total{iter,1}{i,ii}(t+1)./bid_vol_plus_num{iter,1}(i,ii),'r*','MarkerSize',10)
% % %         hold on
% % %         plot(0,ask_vol_plus_total{iter,1}{i,ii}(t+1)./ask_vol_plus_num{iter,1}(i,ii),'b*','MarkerSize',10)
% % %         hold on
% % % %         plot([-t:-1,1:t],bid_vol_plus_total{iter,1}{i,ii}([1:t,t+2:2*t+1])./bid_vol_plus_num{iter,1}(i,ii),'ro')
% % % %         hold on
% % % %         plot([-t:-1,1:t],ask_vol_plus_total{iter,1}{i,ii}([1:t,t+2:2*t+1])./ask_vol_plus_num{iter,1}(i,ii),'bo')
% % %     %     legend('best bid depth','best ask depth')
% % %         xlabel('lag')
% % %         ylabel('Volume')
% % %         if iii==1
% % %             title('Type 1')
% % %         elseif iii==2
% % %             title('Type 7')
% % %         elseif iii==3
% % %             title('Type 11')
% % %         end
% % % %         if i==1
% % % %             title('Type A, minus initiated market order')
% % % %         elseif i==2
% % % %             title('Type B, minus initiated market order')
% % % %         elseif i==3
% % % %             title('Type C, minus initiated market order')
% % % %         elseif i==4
% % % %             title('Type D, minus initiated market order')
% % % %         elseif i==5
% % % %             title('Type E, minus initiated market order')
% % % %         elseif i==6
% % % %             title('Type F, minus initiated market order')
% % % %         elseif i==7
% % % %             title('Type G, minus initiated market order')
% % % %         end
% % %     end
% % % %     saveas(f1,sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/f1/%d_%d_%d_%d_depth_plusside.png',firm,iter,id,ii))
% % % end
% % % 
% % % 
% % % 
% % % %% order book gap
% % % iter=iter_set;
% % % % ii=1;
% % % for ii=1:1
% % %     f1=figure;
% % %     set(gcf,'color','w')
% % % %     set(f1,'Visible','off')
% % % %     set(f1, 'PaperUnits', 'inches');
% % % %     x_width=14 ;y_width=7;
% % % %     set(f1, 'PaperPosition', [0 0 x_width y_width]);
% % %     for i=1:16
% % %         subplot(4,4,i)
% % %         plot(-t:t,bid_gap_plus_total{iter,1}{i,ii}./bid_gap_plus_num{iter,1}(i,ii),'r.','LineWidth',2)
% % %         hold on
% % %         plot(-t:t,ask_gap_plus_total{iter,1}{i,ii}./ask_gap_plus_num{iter,1}(i,ii),'b.','LineWidth',2)
% % %         hold on
% % %         plot(0,bid_gap_plus_total{iter,1}{i,ii}(t+1)./bid_gap_plus_num{iter,1}(i,ii),'r*','MarkerSize',10)
% % %         hold on
% % %         plot(0,ask_gap_plus_total{iter,1}{i,ii}(t+1)./ask_gap_plus_num{iter,1}(i,ii),'b*','MarkerSize',10)
% % % %         hold on
% % % %         plot([-t:-1,1:t],bid_gap_plus_total{iter,1}{i,ii}([1:t,t+2:2*t+1])./bid_gap_plus_num{iter,1}(i,ii),'ro')
% % % %         hold on
% % % %         plot([-t:-1,1:t],ask_gap_plus_total{iter,1}{i,ii}([1:t,t+2:2*t+1])./ask_gap_plus_num{iter,1}(i,ii),'bo')
% % %     %     legend('best bid gap','best ask gap')
% % %         xlabel('lag')
% % %         ylabel('Gap')
% % % %         if i==1
% % % %             title('Type A, minus initiated market order')
% % % %         elseif i==2
% % % %             title('Type B, minus initiated market order')
% % % %         elseif i==3
% % % %             title('Type C, minus initiated market order')
% % % %         elseif i==4
% % % %             title('Type D, minus initiated market order')
% % % %         elseif i==5
% % % %             title('Type E, minus initiated market order')
% % % %         elseif i==6
% % % %             title('Type F, minus initiated market order')
% % % %         elseif i==7
% % % %             title('Type G, minus initiated market order')
% % % %         end
% % %     end
% % % %     saveas(f1,sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/f1/%d_%d_%d_%d_gap_plusside.png',firm,iter,id,ii))
% % % 
% % %     f1=figure;
% % %     set(gcf,'color','w')
% % % %     set(f1,'Visible','off')
% % %     set(f1, 'PaperUnits', 'inches');
% % %     x_width=14 ;y_width=7;
% % %     set(f1, 'PaperPosition', [0 0 x_width y_width]);
% % %     for i=1:NI
% % %         subplot(4,5,i)
% % %         plot(-t:t,bid_gap_minus_total{iter,1}{i,ii}./bid_gap_minus_num{iter,1}(i,ii),'r','LineWidth',2)
% % %         hold on
% % %         plot(-t:t,ask_gap_minus_total{iter,1}{i,ii}./ask_gap_minus_num{iter,1}(i,ii),'b','LineWidth',2)
% % %         hold on
% % %         plot(0,bid_gap_minus_total{iter,1}{i,ii}(t+1)./bid_gap_minus_num{iter,1}(i,ii),'r*','MarkerSize',10)
% % %         hold on
% % %         plot(0,ask_gap_minus_total{iter,1}{i,ii}(t+1)./ask_gap_minus_num{iter,1}(i,ii),'b*','MarkerSize',10)
% % %         hold on
% % %         plot([-t:-1,1:t],bid_gap_minus_total{iter,1}{i,ii}([1:t,t+2:2*t+1])./bid_gap_minus_num{iter,1}(i,ii),'ro')
% % %         hold on
% % %         plot([-t:-1,1:t],ask_gap_minus_total{iter,1}{i,ii}([1:t,t+2:2*t+1])./ask_gap_minus_num{iter,1}(i,ii),'bo')
% % %     %     legend('best bid gap','best ask gap')
% % %         xlabel('lag')
% % %         ylabel('Gap')
% % %         if i==1
% % %             title('Type A, plus initiated market order')
% % %         elseif i==2
% % %             title('Type B, plus initiated market order')
% % %         elseif i==3
% % %             title('Type C, plus initiated market order')
% % %         elseif i==4
% % %             title('Type D, plus initiated market order')
% % %         elseif i==5
% % %             title('Type E, plus initiated market order')
% % %         elseif i==6
% % %             title('Type F, plus initiated market order')
% % %         elseif i==7
% % %             title('Type G, plus initiated market order')
% % %         end
% % %     end
% % % %     saveas(f1,sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/f1/%d_%d_%d_%d_gap_minusside.png',firm,iter,id,ii))
% % % end
% % % 
% % % 
% % % 
% % % %% virtual shock (1 figure)
% % % iter=iter_set;
% % % f1=figure;
% % % set(gcf,'color','w')
% % % % set(f1,'Visible','off')
% % % % set(f1, 'PaperUnits', 'inches');
% % % % x_width=14 ;y_width=7;
% % % % set(f1, 'PaperPosition', [0 0 x_width y_width]);
% % % 
% % % ii=3;
% % % vi = [1000, 5000, 10000, 20000];
% % % % refer(1) = bid_price_plus_total{iter,1}{1,1}(1)./bid_price_plus_num{iter,1}(1,1);
% % % 
% % % 
% % % for i=1:NI
% % %     if i==1
% % %         ind = 1;
% % %     elseif i>1 && i<4
% % %         ind = 2;
% % %         col = (i-1)/3;
% % %     elseif i>3 && i<7
% % %         ind = 3;
% % %         col = (i-3)/4;
% % %     elseif i>6 && i<11
% % %         ind = 4;
% % %         col = (i-6)/5;
% % %     elseif i>10
% % %         ind = 5;
% % %         col = (i-10)/6;
% % %     end
% % %     
% % %     if i==1 || i==2 || i==4 || i==7 || i==11
% % %         flag_g = 1;
% % %     else
% % %         flag_g = 2;
% % %     end
% % %     
% % %     subplot(1,5,ind)
% % % 
% % %     if flag_g==1
% % %         plot(-t:t,bid_vol_price_plus_avg_ta{iter,1}{i,ii},'r-','LineWidth',1)
% % %         hold on
% % %         plot(-t:t,ask_vol_price_plus_avg_ta{iter,1}{i,ii},'b-','LineWidth',1)
% % %         hold on
% % %         plot(0,bid_vol_price_plus_avg_ta{iter,1}{i,ii}(t+1),'r*','MarkerSize',10)
% % %         hold on
% % %         plot(0,ask_vol_price_plus_avg_ta{iter,1}{i,ii}(t+1),'b*','MarkerSize',10)
% % %     elseif flag_g==2
% % %         plot(-t:t,bid_vol_price_plus_avg_ta{iter,1}{i,ii},'-','Color',[1 col col],'LineWidth',1)
% % %         hold on
% % %         plot(-t:t,ask_vol_price_plus_avg_ta{iter,1}{i,ii},'-','Color',[col col 1],'LineWidth',1)
% % %         hold on
% % %         plot(0,bid_vol_price_plus_avg_ta{iter,1}{i,ii}(t+1),'o','Color',[1 col col],'MarkerSize',10)
% % %         hold on
% % %         plot(0,ask_vol_price_plus_avg_ta{iter,1}{i,ii}(t+1),'o','Color',[col col 1],'MarkerSize',10)
% % %     end
% % % %         set(gca,'xscale','log')
% % % %         set(gca,'yscale','log')
% % % %         hold on
% % % %     hold on
% % % %     plot([-t:1:-1,1:t],spread_minus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_minus_num{iter,1}(i,1),'b.')
% % % %     legend('minus initiated market order','plus initiated market order')
% % %     xlabel('lag')
% % %     ylabel('Price Shock')
% % % %         ylim([0.6 18])
% % % end
% % % 
% % % % iter=iter_set;
% % % % f2=figure;
% % % % set(gcf,'color','w')
% % % % % set(f1,'Visible','off')
% % % % % set(f1, 'PaperUnits', 'inches');
% % % % % x_width=14 ;y_width=7;
% % % % % set(f1, 'PaperPosition', [0 0 x_width y_width]);
% % % % 
% % % % ii=3;
% % % % vi = [1000, 5000, 10000, 20000];
% % % % % refer(1) = bid_price_plus_total{iter,1}{1,1}(1)./bid_price_plus_num{iter,1}(1,1);
% % % % 
% % % % 
% % % % for i=1:NI
% % % %     if i>5 && i<10
% % % %         ind = i-4;
% % % %         flag = 1;
% % % %         flag_g = 1;
% % % %     elseif i>9
% % % %         ind = i-8;
% % % %         flag = 1;
% % % %         flag_g = 2;
% % % %     elseif i==1
% % % %         ind = 1;
% % % %         flag = 1;
% % % %         flag_g = 1;
% % % %     elseif i>1 && i<6
% % % %         flag = 0;
% % % %     end
% % % %     
% % % %     if flag==1
% % % %         subplot(1,5,ind)
% % % %         
% % % %         if flag_g==1
% % % %             plot(-t:t,bid_vol_price_minus_avg_ta{iter,1}{i,ii},'r-','LineWidth',1)
% % % %             hold on
% % % %             plot(-t:t,ask_vol_price_minus_avg_ta{iter,1}{i,ii},'b-','LineWidth',1)
% % % %             hold on
% % % %             plot(0,bid_vol_price_minus_avg_ta{iter,1}{i,ii}(t+1),'r*','MarkerSize',10)
% % % %             hold on
% % % %             plot(0,ask_vol_price_minus_avg_ta{iter,1}{i,ii}(t+1),'b*','MarkerSize',10)
% % % %         elseif flag_g==2
% % % %             plot(-t:t,bid_vol_price_minus_avg_ta{iter,1}{i,ii},'r--','LineWidth',1)
% % % %             hold on
% % % %             plot(-t:t,ask_vol_price_minus_avg_ta{iter,1}{i,ii},'b--','LineWidth',1)
% % % %             hold on
% % % %             plot(0,bid_vol_price_minus_avg_ta{iter,1}{i,ii}(t+1),'ro','MarkerSize',10)
% % % %             hold on
% % % %             plot(0,ask_vol_price_minus_avg_ta{iter,1}{i,ii}(t+1),'bo','MarkerSize',10)
% % % %         end
% % % % %         set(gca,'xscale','log')
% % % % %         set(gca,'yscale','log')
% % % % %         hold on
% % % %     %     hold on
% % % %     %     plot([-t:1:-1,1:t],spread_minus_total{iter,1}{i,1}([1:t,t+2:2*t+1])./spread_minus_num{iter,1}(i,1),'b.')
% % % %     %     legend('minus initiated market order','plus initiated market order')
% % % %         xlabel('lag')
% % % %         ylabel('Price Shock')
% % % % %         ylim([3 24])
% % % %     end
% % % % end
% % % 
% % % 
% % % %%
% % % figure100 = figure('Position',[1 1 1100 500],'Color',[1 1 1]);
% % % subplot(2,2,1)
% % % plot(1:10,sin(1:10));grid;
% % % text(-0.1,1.1,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
% % % subplot(2,2,2)
% % % plot(11:20,sin(1:10));grid;
% % % text(-0.1,1.1,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
% % % subplot(2,2,3)
% % % plot(21:30,sin(1:10));grid;
% % % text(-0.1,1.1,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top')
% % % subplot(2,2,4)
% % % plot(31:40,sin(1:10));grid;
% % % text(-0.1,1.1,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top')
% % % AxesHandle=findobj(figure100,'Type','axes');
% % % set(AxesHandle(4),'Position',[0.05,0.5,0.4,0.4]);
% % % set(AxesHandle(3),'Position',[0.5+0.05,0.5,0.4,0.4]);
% % % set(AxesHandle(2),'Position',[0.05,0.05,0.4,0.4]);
% % % set(AxesHandle(1),'Position',[0.5+0.05,0.05,0.4,0.4]);
% % % 
% % % 
% % % 
% % % 
% % % figure100 = figure('Position',[1 1 1100 500],'Color',[1 1 1]);
% % % subplot(1,3,1)
% % % plot(1:10,sin(1:10));grid;
% % % text(-0.1,1.1,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
% % % subplot(1,3,2)
% % % plot(11:20,sin(1:10));grid;
% % % text(-0.1,1.1,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
% % % subplot(1,3,3)
% % % plot(21:30,sin(1:10));grid;
% % % text(-0.1,1.1,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top')
% % % AxesHandle=findobj(figure100,'Type','axes');
% % % set(AxesHandle(3),'Position',[0.05,0.05,0.25,0.85]);
% % % set(AxesHandle(2),'Position',[0.33+0.05,0.05,0.25,0.85]);
% % % set(AxesHandle(1),'Position',[0.66+0.05,0.05,0.25,0.85]);
% % % 
% % % 
