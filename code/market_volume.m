%%
clear;clc;

firm_start = 1;
firm_end = 65;
len_firm = firm_end - firm_start + 1;
iter_set=1; % tick (in order, refer to )
id=1; % spread (1: all, 2: small, 3: large)
day_start=1;
day_end=169;
t=1000;

cd '/home/minyoung/data_minyoung/Research/Price_Impact'
filepath='/home/minyoung/data_minyoung/cm_code/test/data/qth_corr_data/result/LSE_data2';
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/');
load('/home/minyoung/data_minyoung/Research/Price_Impact/date/date.txt')
com_num=load('/home/minyoung/data_minyoung/Research/Price_Impact/top_69.csv');

if day_start<=26 && day_end>=26
    len_day=day_end-day_start+1-1;    
else
    len_day=day_end-day_start+1;
end


tolerance=0.00001; % for comparison of two variable of the double type

iter=0;
NI=12;
de=10;


market_stat = cell(len_day,len_firm);
market_ind_save = cell(len_firm,1);
firm_ind = 1;


for firm=firm_start:firm_end
    firm

    market_ind=cell(NI,2);
    tick_bin_valid_save=[];
    day_ind = 1;
    for day=day_start:day_end
    %     [firm, day]
        market_stat_daily = zeros(NI,2);
    
        if exist(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_AB_sign_%d_%d_%d_%d.mat',com_num(firm),day,id,t),'file') && day~=26 ...
                && exist(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_AB_sign_ind_%d_%d_%d_%d.mat',com_num(firm),day,id,t),'file') ...
                && exist(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_AB_sign_var_%d_%d_%d_%d.mat',com_num(firm),day,id,t),'file') ...
                && exist(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_ABa_sign_%d_%d_%d_%d.mat',com_num(firm),day,id,t),'file') ...
                && exist(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_ABa_sign_ind_%d_%d_%d_%d.mat',com_num(firm),day,id,t),'file') ...
                && exist(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result_over/lmc_flow_event_ver6_AB_over_sign_%d_%d_%d_%d.mat',com_num(firm),day,id,t),'file') ...
                && exist(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result_over/lmc_flow_event_ver6_AB_over_sign_ind_%d_%d_%d_%d.mat',com_num(firm),day,id,t),'file')
            
            load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_AB_sign_%d_%d_%d_%d.mat',com_num(firm),day,id,t),'tick_bin_valid')
            load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_AB_sign_ind_%d_%d_%d_%d.mat',com_num(firm),day,id,t))
            load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_AB_sign_var_%d_%d_%d_%d.mat',com_num(firm),day,id,t),'market')

    %         ABa1 = load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_ABa_sign_%d_%d_%d_%d.mat',com_num(firm),day,id,t));
            ABa2 = load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_ABa_sign_ind_%d_%d_%d_%d.mat',com_num(firm),day,id,t));
            AB_over = load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result_over/lmc_flow_event_ver6_AB_over_sign_ind_%d_%d_%d_%d.mat',com_num(firm),day,id,t));
    %         ABa3 = load(sprintf('/home/minyoung/data_minyoung/Research/Price_Impact/data/qth_corr_data/result/lmc_flow_event_ver6_ABa_sign_var_%d_%d_%d_%d.mat',com_num(firm),day,id,t));  

            for i=1:length(tick_bin_valid)
                ind_minus{10,1}{i,1} = ABa2.ind_minus{1,1}{i,1};
                ind_plus{10,1}{i,1} = ABa2.ind_plus{1,1}{i,1};
                
                ind_minus{11,1}{i,1} = AB_over.ind_minus{1,1}{i,1};
                ind_plus{11,1}{i,1} = AB_over.ind_plus{1,1}{i,1};
                ind_minus{12,1}{i,1} = AB_over.ind_minus{2,1}{i,1};
                ind_plus{12,1}{i,1} = AB_over.ind_plus{2,1}{i,1};
            end
            clear ABa1 ABa2 AB_over
            
            for ii=1:length(tick_bin_valid)
                for i=1:NI
                    market_stat_daily(i,1) = market_stat_daily(i,1) + length(ind_minus{i,1}{ii,1});
                    market_stat_daily(i,2) = market_stat_daily(i,2) + length(ind_plus{i,1}{ii,1});
                end
            end
            
            market_stat{day_ind,firm_ind} = market_stat_daily;
            clear market_stat_daily


            ch=0;
            for ii=1:length(tick_bin_valid_save)
                for i=1:length(tick_bin_valid)
                    if find(tick_bin_valid_save(ii)>tick_bin_valid(i)*(1-tolerance) & tick_bin_valid_save(ii)<tick_bin_valid(i)*(1+tolerance))
                        tick_bin_valid(i)=tick_bin_valid_save(ii);
                        ch=ch+1;
                    end
                end
            end

            tick_bin_valid_save=[tick_bin_valid_save;tick_bin_valid];
            [dummy ii]=unique(tick_bin_valid_save,'first');
            tick_bin_valid_save=tick_bin_valid_save(sort(ii));


            for i=1:NI
                minus_ind=[];
                for ii=1:length(eval(sprintf('ind_minus{%d,1}',i)))
                    minus_ind = [minus_ind; eval(sprintf('ind_minus{%d,1}{%d,1}',i,ii))];
                end
                market_ind{i,1}=[market_ind{i,1};market(minus_ind,:)];

                plus_ind=[];
                for ii=1:length(eval(sprintf('ind_plus{%d,1}',i)))
                    plus_ind = [plus_ind; eval(sprintf('ind_plus{%d,1}{%d,1}',i,ii))];
                end
                market_ind{i,2}=[market_ind{i,2};market(plus_ind,:)];
                clear ind_a ind_b
            end
            clear minus_ind plus_ind

            clear ind_minus ind_plus
            clear ask_gap_minus ask_gap_plus ask_price_minus ask_price_plus ask_vol_minus ask_vol_plus
            clear bid_gap_minus bid_gap_plus bid_price_minus bid_price_plus bid_vol_minus bid_vol_plus spread_minus spread_plus flow_minus flow_plus
            clear bid_vol_price_minus bid_vol_price_plus ask_vol_price_minus ask_vol_price_plus
            
            day_ind = day_ind + 1;
        elseif day==26
        else
            sprintf('Error: empty files: %d, %d',firm,day)
        end
        clear aa ask ask_di ask_gap ask_gap_event ask_price ask_price_event ask_return ask_vol ask_vol_event
        clear bid bid_di bid_gap bid_gap_event bid_price bid_price_event bid_return bid_vol bid_vol_event
        clear cancel flow lc_minus limit market orderflow price price_return spread spread_a spread_b t_time
        clear time_event time_gap flow trade_id update_time_gap_minus update_time_gap_plus
        clear bid_vol_price_minus_avg bid_vol_price_plus_avg ask_vol_price_minus_avg ask_vol_price_plus_avg
        
%         day_ind = day_ind + 1;
    end
    market_ind_save{firm_ind,1} = market_ind;
    firm_ind = firm_ind + 1;
    clear market_ind
end

%%
clear;close all;clc;

% load('market_stat_save.mat')
% load('market_stat_save_70.mat')
load('market_stat_save_65.mat')

%% market_statistics

sum_trading = zeros(len_day,len_firm);

for i=1:len_day
    for ii=1:len_firm
        sum_trading(i,ii) = sum(sum(market_stat{i,ii}));
    end
end

stat_ni = cell(NI,1);
stat_ni_firm = cell(NI,1);
for ni = 1:NI
    for i=1:len_day
        temp = [];
        for ii=1:len_firm
            temp = [temp; market_stat{i,ii}(ni,2)];
            stat_ni_firm{ni,1}(ii,i) = market_stat{i,ii}(ni,2);
        end
        stat_ni{ni,1} = [stat_ni{ni,1}; mean(temp)];
    end
end

for ni = 1:NI
    stat_sum{ni,1} = sum(stat_ni_firm{ni,1}(:,1:30),2);
end

figure
set(gcf,'color','w')
for i=1:NI
    subplot(3,4,i)
    plot(stat_ni{i,1})
    xlim([1 60])
end

figure
set(gcf,'color','w')
a = mean(sum_trading);
a = a'/(8.5*60*60);
bar(1./a)
xlim([0 31])
xlabel('firm')
ylabel('average trading interval (sec)')

%%
% average over firm
market_stat_tot = cell(1,len_firm);
for i=1:len_firm
    market_stat_tot{1,i} = zeros(NI,2);
    for ii=1:len_day
        market_stat_tot{1,i} = market_stat_tot{1,i} + market_stat{ii,i};
    end
    market_stat_tot{1,i} = market_stat_tot{1,i}./(len_day);
end


for i=1:len_firm
    market_stat_firm(i,1) = sum(market_stat_tot{1,i}(10,:));
    market_stat_firm(i,2) = sum(market_stat_tot{1,i}(1,:)+market_stat_tot{1,i}(2,:)+market_stat_tot{1,i}(4,:)+market_stat_tot{1,i}(6,:)+market_stat_tot{1,i}(8,:)+market_stat_tot{1,i}(11,:));
    market_stat_firm(i,3) = sum(market_stat_tot{1,i}(3,:)+market_stat_tot{1,i}(5,:)+market_stat_tot{1,i}(7,:)+market_stat_tot{1,i}(9,:)+market_stat_tot{1,i}(12,:));
end
two_firm = market_stat_firm(:,2:3);
ratio_firm(:,1) = two_firm(:,1)./sum(two_firm,2);
ratio_firm(:,2) = two_firm(:,2)./sum(two_firm,2);

figure
set(gcf,'color','w')
subplot(1,3,1)
bar(ratio_firm*100,'stacked')
ylim([95 100])

ind = 1;
for i=1:len_firm
    market_stat_firm(i,1) = market_stat_tot{1,i}(10,ind);
    market_stat_firm(i,2) = market_stat_tot{1,i}(1,ind)+market_stat_tot{1,i}(2,ind)+market_stat_tot{1,i}(4,ind)+market_stat_tot{1,i}(6,ind)+market_stat_tot{1,i}(8,ind)+market_stat_tot{1,i}(11,ind);
    market_stat_firm(i,3) = market_stat_tot{1,i}(3,ind)+market_stat_tot{1,i}(5,ind)+market_stat_tot{1,i}(7,ind)+market_stat_tot{1,i}(9,ind)+market_stat_tot{1,i}(12,ind);
end
two_firm = market_stat_firm(:,2:3);
ratio_firm(:,1) = two_firm(:,1)./sum(two_firm,2);
ratio_firm(:,2) = two_firm(:,2)./sum(two_firm,2);

subplot(1,3,2)
bar(ratio_firm*100,'stacked')
ylim([95 100])

ind = 2;
for i=1:len_firm
    market_stat_firm(i,1) = market_stat_tot{1,i}(10,ind);
    market_stat_firm(i,2) = market_stat_tot{1,i}(1,ind)+market_stat_tot{1,i}(2,ind)+market_stat_tot{1,i}(4,ind)+market_stat_tot{1,i}(6,ind)+market_stat_tot{1,i}(8,ind)+market_stat_tot{1,i}(11,ind);
    market_stat_firm(i,3) = market_stat_tot{1,i}(3,ind)+market_stat_tot{1,i}(5,ind)+market_stat_tot{1,i}(7,ind)+market_stat_tot{1,i}(9,ind)+market_stat_tot{1,i}(12,ind);
end
two_firm = market_stat_firm(:,2:3);
ratio_firm(:,1) = two_firm(:,1)./sum(two_firm,2);
ratio_firm(:,2) = two_firm(:,2)./sum(two_firm,2);

subplot(1,3,3)
bar(ratio_firm*100,'stacked')
ylim([95 100])

%%
% day
market_stat_tot_day = cell(1,len_day);
for i=1:len_day
    market_stat_tot_day{1,i} = zeros(NI,2);
    for ii=1:len_firm
        market_stat_tot_day{1,i} = market_stat_tot_day{1,i} + market_stat{i,ii};
    end
    market_stat_tot_day{1,i} = market_stat_tot_day{1,i}./(len_firm);
end

ind = 1:2;
clear market_stat_tot_4
for i=1:len_day
    market_stat_day(i,1) = sum(market_stat_tot_day{1,i}(10,ind));
    market_stat_day(i,2) = sum(market_stat_tot_day{1,i}(2,ind)+market_stat_tot_day{1,i}(4,ind)+market_stat_tot_day{1,i}(6,ind)+market_stat_tot_day{1,i}(8,ind)+market_stat_tot_day{1,i}(11,ind));
    market_stat_day(i,3) = sum(market_stat_tot_day{1,i}(3,ind)+market_stat_tot_day{1,i}(5,ind)+market_stat_tot_day{1,i}(7,ind)+market_stat_tot_day{1,i}(9,ind)+market_stat_tot_day{1,i}(12,ind));
    market_stat_day(i,4) = sum(market_stat_tot_day{1,i}(1,ind));
end


two_day = market_stat_day(:,2:3);
ratio_day = two_day;
ratio_day(:,1) = two_day(:,1)./sum(two_day,2);
ratio_day(:,2) = two_day(:,2)./sum(two_day,2);


addpath('/home/minyoung/data_minyoung/Research/Matlab_function')

market_stat_day_new(1:25,:) = market_stat_day(1:25,:);
market_stat_day_new(26,:) = [0, 0, 0, 0];
market_stat_day_new(27:169,:) = market_stat_day(26:end,:);

market_stat_mv = moving_sum(market_stat_day_new,10,1,1);


addpath('/home/minyoung/data_minyoung/Research/Matlab_function/linspecer')
C2 = linspecer(7);

figure1 = figure(1);
set(gcf,'color','w')
plot(market_stat_mv(:,1),'-','color',[0 0 0],'linewidth',2)
hold on
plot(market_stat_mv(:,4),':','color',C2(2,:),'linewidth',2)
hold on
plot(market_stat_mv(:,2),'--','color',C2(4,:),'linewidth',2)
hold on
plot(market_stat_mv(:,3),'-.','color',C2(1,:),'linewidth',2)

hold on
plot([22 22],[10^2 4*10^4],'k--')
set(gca,'yscal','log')
legend('Type Zero','Type One','Type A','Type B')
xlabel('Date','fontsize',18)
ylabel('Number of market order','fontsize',18)
ylim([10^2 4*10^4])

date_seq = load('./date/date.txt');
date_seq = date_seq(10:end);
date_str = cell(160,1);
for i=1:160
    date_str{i,1} = num2str(date_seq(i,1));
end
date_first = [20080815;20080915;20081015;20081114;20081215;20090115;20090216;20090316];
for i=1:length(date_first)
    date_first_ind(i,1) = find(date_seq == date_first(i));
end
date_str = date_str';
date_str_tick = date_str(date_first_ind);

set(gca,'xtick',date_first_ind)
set(gca,'xticklabel',date_str_tick)

xtickangle(35)



figure
set(gcf,'color','w')
plotyy(1:size(two_day,1),two_day(:,1),1:size(two_day,1),two_day(:,2))

figure
set(gcf,'color','w')
bar(ratio_day*100,'stacked')
ylim([90 100])
%% ccdf of market_ind for all firms
addpath('/home/minyoung/data_minyoung/Research/Matlab_function');
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/linspecer')

C2 = linspecer(10);
binsz = 100;

firm = 1; % AZN
market_ind = market_ind_save{firm,1};
figure1 = figure;
set(gcf,'color','w')
ind(:,1) = [10;1;2;4;6;8];
ind(:,2) = [10;3;5;7;9;10];
subplot(1,3,1)
for i=1:6
    [a b c d] = pcdf(market_ind{ind(i,1),2}(:,2),binsz);

    if i==1
        plot(a,d,'--','color',[0 0 0])
        hold on
    else
        plot(a,d,'color',C2(i-1,:))
        hold on
    end
    if i==6
        set(gca,'xscale','log')
        set(gca,'yscale','log')
    end
end
xlabel('V')
ylabel('CCDF')
xlim([10 inf])
text(-0.15,1.17,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
legend('\Delta a_0=0 (type Zero)','\Delta a_0=1 (type A)','\Delta a_0=2 (type A)','\Delta a_0=3 (type A)','\Delta a_0=4 (type A)','\Delta a_0=5 (type A)')
legend boxoff

subplot(1,3,2)
for i=1:5
    [a b c d] = pcdf(market_ind{ind(i,2),2}(:,2),binsz);
    
    if i==1
        plot(a,d,'--','color',[0 0 0])
        hold on
    else
        plot(a,d,'color',C2(i+4,:))
        hold on
    end
    if i==5
        set(gca,'xscale','log')
        set(gca,'yscale','log')
    end
end
xlabel('V')
ylabel('CCDF')
xlim([10 inf])
text(-0.15,1.17,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
legend('\Delta a_0=0 (type Zero)','\Delta a_0=2 (type B)','\Delta a_0=3 (type B)','\Delta a_0=4 (type B)','\Delta a_0=5 (type B)')
legend boxoff

ind_t(:,1) = [10;1;2;4;6;8;3;5;7;9];
subplot(1,3,3)
for i=1:10
    [a b c d] = pcdf(market_ind{i,2}(:,2)/std(market_ind{i,2}(:,2)),binsz);
    
    if i==1
        plot(a,d,'--','color',[0 0 0])
        hold on
    else
        plot(a,d,'color',C2(i-1,:))
        hold on
    end
    if i==10
        set(gca,'xscale','log')
        set(gca,'yscale','log') 
    end
end
xlabel('V / \sigma')
ylabel('CCDF')
xlim([(10^-1)/5 inf])
text(-0.15,1.17,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top')
% legend('\Delta a_0=0','\Delta a_0=1','\Delta a_0=2','\Delta a_0=3','\Delta a_0=4',...
%     '\Delta a_0=5','\Delta a_0=2','\Delta a_0=3','\Delta a_0=4','\Delta a_0=5')
% legend boxoff

dx=0.05;
dy=0.05;
x=(1-5*dx)/2.9;
y=(1-5*dy)/1.1;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure1,'Type','axes');
set(AxesHandle(3),'Position',[dx+dx/10,dy+dy/2+1.8*dy,x,y]);
set(AxesHandle(2),'Position',[x+2.5*dx+dx/10,dy+dy/2+1.8*dy,x,y]);
set(AxesHandle(1),'Position',[2*x+4*dx+dx/10,dy+dy/2+1.8*dy,x,y]);
set(gcf,'Position',[100, 100, 1000, 330])


%% ccdf of market_ind for all firms  18.12.14 (minus + plus)

addpath('/home/minyoung/data_minyoung/Research/Matlab_function');
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/linspecer')

C2 = linspecer(10);
binsz = 100;

firm = 1; % 1:HSBC, 12:AZN
market_ind = market_ind_save{firm,1};
figure1 = figure;
set(gcf,'color','w')
ind(:,1) = [10;1;2;4;6;8;11];
ind(:,2) = [10;3;5;7;9;12;10]; % 10 unused
subplot(1,2,1)
for i=1:7
    if i == 1
        [a b c d] = pcdf([market_ind{ind(i,1),2}(:,2); market_ind{ind(i,1),1}(:,1)],binsz);
    elseif i == 2
        [a b c d] = pcdf([market_ind{ind(i,1),2}(:,2);market_ind{ind(i,1),1}(:,1)],binsz);
    else
        [a b c d] = pcdf([market_ind{ind(i,1),2}(:,2); market_ind{ind(i-1,2),2}(:,2);market_ind{ind(i,1),1}(:,1); market_ind{ind(i-1,2),1}(:,1)],binsz);
    end

    if i==1
        plot(a,d,'--','color',[0 0 0])
        hold on
    elseif i == 2
        plot(a,d,':','color',[2/10 2/10 2/10])
        hold on   
    else
%         plot(a,d,'color',C2(i-1,:))
        axis = plot(a,d,'color',[i/10 i/10 i/10]);
        hold on
    end
    if i==7
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        set(gca,'Layer','top')
    end
end
xlabel('V')
ylabel('CCDF')
xlim([5*10^2 inf])
text(-0.15,1.17,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
legend('d_a(or d_b)=0  (Type Zero)','d_a(or d_b)=1  (Type One)','d_a(or d_b)=2  (Type A + B)','d_a(or d_b)=3  (Type A + B)','d_a(or d_b)=4  (Type A + B)','d_a(or d_b)=5  (Type A + B)','d_a(or d_b)>5  (Type A + B)')
legend boxoff

subplot(1,2,2)
for i=1:7
    %[a b c d] = pcdf(market_ind{ind(i,1),2}(:,2),binsz); % 19.04.25 check
    [a b c d] = pcdf([market_ind{ind(i,1),2}(:,2); market_ind{ind(i,1),1}(:,1)],binsz);
    
    if i==1
        plot(a,d,'--','color',[0 0 0])
        hold on
    elseif i == 2
        plot(a,d,':','color',[2/10 2/10 2/10])
        hold on
    else
%         plot(a,d,'color',C2(i-1,:))
        plot(a,d,'color',[i/10 i/10 i/10])
%         plot(a,d)
        hold on
    end
    if i==6
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        set(gca,'Layer','top')
    end
end
xlabel('V')
ylabel('CCDF')
xlim([5*10^2 inf])
text(-0.15,1.17,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
legend('d_a(or d_b)=0  (Type Zero)','d_a(or d_b)=1  (Type One)','d_a(or d_b)=2  (Type A)','d_a(or d_b)=3  (Type A)','d_a(or d_b)=4  (Type A)','d_a(or d_b)=5  (Type A)','d_a(or d_b)>5  (Type A)')
legend boxoff

dx=0.05;
dy=0.05;
x=(1-5*dx)/1.9;
y=(1-5*dy)/1.1;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure1,'Type','axes');
set(AxesHandle(2),'Position',[dx+dx/2,dy+dy/2+1.8*dy,x,y]);
set(AxesHandle(1),'Position',[x+2.5*dx+dx,dy+dy/2+1.8*dy,x,y]);
set(gcf,'Position',[100, 100, 850, 400])


%% table of market_ind for all firms  19.4.25 (minus + plus)

addpath('/home/minyoung/data_minyoung/Research/Matlab_function');
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/linspecer')

C2 = linspecer(10);
binsz = 100;

firm = 1; % 1:HSBC, 12:AZN

ind(:,1) = [10;1;2;4;6;8;11];
ind(:,2) = [10;3;5;7;9;12;10]; % 10 unused

% consider all level of price change
% vol_ratio = cell(70,1);
% 
% for firm = 1:70
%     
%     market_ind = market_ind_save{firm,1};
%     
%     vol_ratio{firm,1} = zeros(12,1);
%     
%     for i=1:5
%         temp = [market_ind{ind(i+1,2),2}(:,2); market_ind{ind(i+1,2),1}(:,1)];
%         vol_ratio{firm,1}(i,1) = length(find(temp>10000))/length(temp);
%     end
%     
% 
%     for i=1:7
%         temp = [market_ind{ind(i,1),2}(:,2); market_ind{ind(i,1),1}(:,1)];
%         vol_ratio{firm,1}(i+5,1) = length(find(temp>10000))/length(temp);
%     end
% end

% consider only type A and B, and Zero
vol_ratio = zeros(70,3);
for firm = 1:70
    
    market_ind = market_ind_save{firm,1};
    
    temp = market_ind{ind(1,1),2}(:,2); market_ind{ind(1,1),1}(:,1);
    vol_ratio(firm,1) = length(find(temp>10000))/length(temp);
    
    temp = [];
    for i=1:6
        temp = [temp; market_ind{ind(i+1,1),2}(:,2); market_ind{ind(i+1,1),1}(:,1)];
    end
    vol_ratio(firm,2) = length(find(temp>10000))/length(temp);
    
    temp = [];
    for i=1:5
        temp = [temp; market_ind{ind(i+1,2),2}(:,2); market_ind{ind(i+1,2),1}(:,1)];
    end
    vol_ratio(firm,3) = length(find(temp>10000))/length(temp);
end


%% ccdf of market_ind for all firms  18.12.15 (minus + plus, cumulative)

addpath('/home/minyoung/data_minyoung/Research/Matlab_function');
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/linspecer')

C2 = linspecer(10);
binsz = 100;

firm = 1; % AZN
market_ind = market_ind_save{firm,1};
figure1 = figure;
set(gcf,'color','w')
ind(:,1) = [10;1;2;4;6;8;11];
ind(:,2) = [10;3;5;7;9;12;10]; % 10 unused
subplot(1,2,1)
for i=1:7
    if i == 1
        [a b c d] = pcdf([market_ind{ind(i,1),2}(:,2); market_ind{ind(i,1),1}(:,1)],binsz);
    elseif i == 2
        [a b c d] = pcdf([market_ind{ind(i,1),2}(:,2);market_ind{ind(i,1),1}(:,1)],binsz);
    else
        [a b c d] = pcdf([market_ind{ind(i,1),2}(:,2); market_ind{ind(i-1,2),2}(:,2);market_ind{ind(i,1),1}(:,1); market_ind{ind(i-1,2),1}(:,1)],binsz);
    end

    if i==1
        plot(a,d,'--','color',[0 0 0])
        hold on
    else
%         plot(a,d,'color',C2(i-1,:))
        plot(a,d,'color',[i/10 i/10 i/10])
        hold on
    end
    if i==7
        set(gca,'xscale','log')
        set(gca,'yscale','log')
    end
end
xlabel('V')
ylabel('CCDF')
xlim([10^2 inf])
text(-0.15,1.17,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
legend('d_a(or d_b)=0 (Type Zero)','d_a(or d_b)=1 (Type A)','d_a(or d_b)=2 (Type A + B)','d_a(or d_b)=3 (Type A + B)','d_a(or d_b)=4 (Type A + B)','d_a(or d_b)=5 (Type A + B)','d_a(or d_b)>5 (Type A + B)')
legend boxoff

subplot(1,2,2)
for i=1:7
    [a b c d] = pcdf(market_ind{ind(i,1),2}(:,2),binsz);
    
    if i==1
        plot(a,d,'--','color',[0 0 0])
        hold on
    else
%         plot(a,d,'color',C2(i-1,:))
        plot(a,d,'color',[i/10 i/10 i/10])
%         plot(a,d)
        hold on
    end
    if i==6
        set(gca,'xscale','log')
        set(gca,'yscale','log')
    end
end
xlabel('V')
ylabel('CCDF')
xlim([10^2 inf])
text(-0.15,1.17,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
legend('d_a(or d_b)=0 (Type Zero)','d_a(or d_b)=1 (Type A)','d_a(or d_b)=2 (Type A)','d_a(or d_b)=3 (Type A)','d_a(or d_b)=4 (Type A)','d_a(or d_b)=5 (Type A)','d_a(or d_b)>5 (Type A)')
legend boxoff

dx=0.05;
dy=0.05;
x=(1-5*dx)/1.9;
y=(1-5*dy)/1.1;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure1,'Type','axes');
set(AxesHandle(2),'Position',[dx+dx/2,dy+dy/2+1.8*dy,x,y]);
set(AxesHandle(1),'Position',[x+2.5*dx+dx,dy+dy/2+1.8*dy,x,y]);
set(gcf,'Position',[100, 100, 800, 350])

%%

market_vol = zeros(1,12);
for firm = 1:1
    for i = 1:12
        market_vol(firm,i) = length([market_ind_save{firm,1}{i,1}(:,1); market_ind_save{firm,1}{i,2}(:,2)]); 
    end
end

market_vol_mean = mean(market_vol,1);
market_vol_mean_per = 100*market_vol_mean/sum(market_vol_mean);

per_ZAB = [sum(market_vol_mean_per(10)); sum(market_vol_mean_per([1,2,4,6,8,11])); sum(market_vol_mean_per([3,5,7,9,12]))];

market_vol_price = zeros(1,7);
market_vol_price(:,1) = market_vol(:,10);
market_vol_price(:,2) = market_vol(:,1);
market_vol_price(:,3) = sum(market_vol(:,[2,3]),2);
market_vol_price(:,4) = sum(market_vol(:,[4,5]),2);
market_vol_price(:,5) = sum(market_vol(:,[6,7]),2);
market_vol_price(:,6) = sum(market_vol(:,[8,9]),2);
market_vol_price(:,7) = sum(market_vol(:,[11,12]),2);

per_price_pz = 100 * market_vol_price ./ sum(market_vol_price,2);
per_price_pz_mean = mean(per_price_pz,1);
per_price_pz_mean_cumsum = cumsum(per_price_pz_mean);

per_price_ez = 100 * market_vol_price(:,2:end) ./ sum(market_vol_price(:,2:end),2);
per_price_ez_mean = mean(per_price_ez,1);
per_price_ez_std = std(per_price_ez,1);
per_price_ez_mean_cumsum = cumsum(mean(per_price_ez,1));


%% fitting for normalized market order volume
addpath('/home/minyoung/data_minyoung/Research/Matlab_function');

binsz = 100;
alpha = cell(70,1);

for firm = 1:70
%     firm = 1; % AZN
    market_ind = market_ind_save{firm,1};
    alpha{firm,1} = zeros(10,1);

    for i=1:10
        [firm,i]
        [alpha{firm,1}(i,1),temp1,temp2] = plfit(market_ind{i,2}(:,2)/std(market_ind{i,2}(:,2)));
        clear temp1 temp2
    end
end

alphas = [];

for i=1:70
    for ii=1:10
        if alpha{i,1}(ii,1) ~= 0
            alphas = [alphas; alpha{i,1}(ii,1)];
        end
    end
end


%%

figure
set(gcf,'color','w')
for firm=1:1 %len_firm
    firm
    market_ind = market_ind_save{firm,1};
    
    figure
    set(gcf,'color','w')

    a00_down = market_ind{10,1}(:,1);
    a00_up = market_ind{10,2}(:,2);
    a1_down = [market_ind{1,1}(:,1);market_ind{2,1}(:,1);market_ind{4,1}(:,1);market_ind{6,1}(:,1);market_ind{8,1}(:,1)];
    a1_up = [market_ind{1,2}(:,2);market_ind{2,2}(:,2);market_ind{4,2}(:,2);market_ind{6,2}(:,2);market_ind{8,2}(:,2)];
    a2_down = [market_ind{3,1}(:,1);market_ind{5,1}(:,1);market_ind{7,1}(:,1);market_ind{9,1}(:,1)];
    a2_up = [market_ind{3,2}(:,2);market_ind{5,2}(:,2);market_ind{7,2}(:,2);market_ind{9,2}(:,2)];
    binsz = 100;

    [a b c d] = pcdf(a00_up,binsz);
    subplot(1,2,1)
    scatter(a,d,30,[0 0 0],'^')
    hold on

    [a b c d] = pcdf(a00_down,binsz);
    subplot(1,2,1)
    scatter(a,d,30,[0 0 0],'v')
    hold on

    [a b c d] = pcdf(a1_up,binsz);
    subplot(1,2,1)
    scatter(a,d,30,[1 0 0],'filled')
    hold on

    [a b c d] = pcdf(a1_down,binsz);
    subplot(1,2,1)
    scatter(a,d,30,[1 0 0])
    hold on

    [a b c d] = pcdf(a2_up,binsz);
    subplot(1,2,1)
    scatter(a,d,30,[0 0 1],'filled')
    hold on

    [a b c d] = pcdf(a2_down,binsz);
    subplot(1,2,1)
    scatter(a,d,30,[0 0 1])
    hold on
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylabel('CCDF, All')
    xlabel('trading volume')
    ylim([5*10^-6 1])

    [a b c d] = pcdf(a00_up/mean(a00_up),binsz);
    subplot(1,2,2)
    scatter(a,d,30,[0 0 0],'^')
    hold on

    [a b c d] = pcdf(a00_down/mean(a00_up),binsz);
    subplot(1,2,2)
    scatter(a,d,30,[0 0 0],'v')
    hold on

    [a b c d] = pcdf(a1_up/mean(a1_up),binsz);
    subplot(1,2,2)
    scatter(a,d,30,[1 0 0],'filled')
    hold on

    [a b c d] = pcdf(a1_down/mean(a1_down),binsz);
    subplot(1,2,2)
    scatter(a,d,30,[1 0 0])
    hold on

    [a b c d] = pcdf(a2_up/mean(a2_up),binsz);
    subplot(1,2,2)
    scatter(a,d,30,[0 0 1],'filled')
    hold on

    [a b c d] = pcdf(a2_down/mean(a2_down),binsz);
    subplot(1,2,2)
    scatter(a,d,30,[0 0 1])
    hold on

    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylabel('CCDF, All')
    xlabel('normalized trading volume')
    ylim([5*10^-6 1])
    
    clear market_ind
end


