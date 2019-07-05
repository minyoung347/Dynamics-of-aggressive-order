%%
clear;clc;

NI = 4;
win_size = 1;
mv_size = 1;
burst_ind = 1;

N = 65;

addpath('/home/minyoung/data_minyoung/Research/Matlab_function')
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/linspecer')
load(sprintf('./burst_cal/burst_individual_1906_%d_%d_%d.mat',win_size,mv_size,burst_ind))
% load(sprintf('./burst_cal/burst_result_1906_%d_%d_%d.mat',win_size,mv_size,burst_ind))


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
                
                ci_burst_minus(i,ii) = 0;
                ci_mem_minus(i,ii) = 0;

                ci_burst_plus(i,ii) = 0;
                ci_mem_plus(i,ii) = 0;
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
            
            ci_burst_minus(i,ii) = tinv(0.975, N-1) * std(burst{i,1}(:,ii)) / sqrt(N);
            ci_mem_minus(i,ii) = tinv(0.975, N-1) * std(mem_coef{i,1}(:,ii)) / sqrt(N);

            ci_burst_plus(i,ii) = tinv(0.975, N-1) * std(burst{i,2}(:,ii)) / sqrt(N);
            ci_mem_plus(i,ii) = tinv(0.975, N-1) * std(mem_coef{i,2}(:,ii)) / sqrt(N);
        end
    end
end

%% statistics for 30 firms (using averaged days)


addpath('/home/minyoung/data_minyoung/Research/Matlab_function/boxplot')
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/sigstar')

alpha_ind = 0; % ks-test, 0 for alpha = 0.01, 1 for alpha = 0.05

len = 3;

% burst_b = zeros(30,3);
% burst_c = zeros(30,3);
% burst_a = zeros(30,3);

burst_b = [];
burst_c = [];
burst_a = [];

% for i=1:len
%     burst_b = burst_b + burst{i,1}/len;
%     burst_c = burst_c + burst{i+30,1}/len;
%     burst_a = burst_a + burst{i+169-win_size+1-len,1}/len;
% end
% burst_b = (burst{1,1} + burst{2,1} + burst{3,1})/3;
% burst_c = (burst{4,1} + burst{5,1} + burst{6,1})/3;
% burst_a = (burst{14,1} + burst{15,1} + burst{16,1})/3;

burst_b = [burst{1,1}; burst{2,1}; burst{3,1}];
burst_c = [burst{4,1}; burst{5,1}; burst{6,1}];
burst_a = [burst{14,1}; burst{15,1}; burst{16,1}];

% burst_b = burst{1,1};
% burst_c = burst{2,1};
% burst_a = burst{5,1};


burst_type_z(:,1) = burst_b(:,1);
burst_type_z(:,2) = burst_c(:,1);
burst_type_z(:,3) = burst_a(:,1);

burst_type_a(:,1) = burst_b(:,2);
burst_type_a(:,2) = burst_c(:,2);
burst_type_a(:,3) = burst_a(:,2);

burst_type_b(:,1) = burst_b(:,3);
burst_type_b(:,2) = burst_c(:,3);
burst_type_b(:,3) = burst_a(:,3);




test_b(1,1) = kstest2(burst_b(:,1),burst_b(:,2),'Alpha',0.01);
test_b(2,1) = kstest2(burst_b(:,1),burst_b(:,2),'Alpha',0.05);
test_b(3,1) = kstest2(burst_b(:,1),burst_b(:,3),'Alpha',0.01);
test_b(4,1) = kstest2(burst_b(:,1),burst_b(:,3),'Alpha',0.05);
test_b(5,1) = kstest2(burst_b(:,2),burst_b(:,3),'Alpha',0.01);
test_b(6,1) = kstest2(burst_b(:,2),burst_b(:,3),'Alpha',0.05);

test_c(1,1) = kstest2(burst_c(:,1),burst_c(:,2),'Alpha',0.01);
test_c(2,1) = kstest2(burst_c(:,1),burst_c(:,2),'Alpha',0.05);
test_c(3,1) = kstest2(burst_c(:,1),burst_c(:,3),'Alpha',0.01);
test_c(4,1) = kstest2(burst_c(:,1),burst_c(:,3),'Alpha',0.05);
test_c(5,1) = kstest2(burst_c(:,2),burst_c(:,3),'Alpha',0.01);
test_c(6,1) = kstest2(burst_c(:,2),burst_c(:,3),'Alpha',0.05);

test_a(1,1) = kstest2(burst_a(:,1),burst_a(:,2),'Alpha',0.01);
test_a(2,1) = kstest2(burst_a(:,1),burst_a(:,2),'Alpha',0.05);
test_a(3,1) = kstest2(burst_a(:,1),burst_a(:,3),'Alpha',0.01);
test_a(4,1) = kstest2(burst_a(:,1),burst_a(:,3),'Alpha',0.05);
test_a(5,1) = kstest2(burst_a(:,2),burst_a(:,3),'Alpha',0.01);
test_a(6,1) = kstest2(burst_a(:,2),burst_a(:,3),'Alpha',0.05);


% figure
% subplot(1,3,1)
% bplot(burst_b)
% set(gca,'XTick',1:3,'XTickLabel',{'Zero','A','B'})
% H = sigstar({[1,2], [1,3], [2,3]},[0.01, 0.01, 0.01]);
% set(H,'color','k')
% ylim([.1 .55])
% 
% subplot(1,3,2)
% bplot(burst_c)
% set(gca,'XTick',1:3,'XTickLabel',{'Zero','A','B'})
% H = sigstar({[1,2], [1,3], [2,3]},[0.01, 0.01, nan]);
% set(H,'color','k')
% ylim([.1 .55])
% 
% subplot(1,3,3)
% bplot(burst_a)
% set(gca,'XTick',1:3,'XTickLabel',{'Zero','A','B'})
% H = sigstar({[1,2], [1,3], [2,3]},[0.01, 0.01, 0.01]);
% set(H,'color','k')
% ylim([.1 .55])

test_type_z = zeros(6,2);
test_type_a = zeros(6,2);
test_type_b = zeros(6,2);

[test_type_z(1,1) test_type_z(1,2)] = kstest2(burst_type_z(:,1),burst_type_z(:,2),'Alpha',0.01);
[test_type_z(2,1) test_type_z(2,2)] = kstest2(burst_type_z(:,1),burst_type_z(:,2),'Alpha',0.05);
[test_type_z(3,1) test_type_z(3,2)] = kstest2(burst_type_z(:,1),burst_type_z(:,3),'Alpha',0.01);
[test_type_z(4,1) test_type_z(4,2)] = kstest2(burst_type_z(:,1),burst_type_z(:,3),'Alpha',0.05);
[test_type_z(5,1) test_type_z(5,2)] = kstest2(burst_type_z(:,2),burst_type_z(:,3),'Alpha',0.01);
[test_type_z(6,1) test_type_z(6,2)] = kstest2(burst_type_z(:,2),burst_type_z(:,3),'Alpha',0.05);

[test_type_a(1,1) test_type_a(1,2)] = kstest2(burst_type_a(:,1),burst_type_a(:,2),'Alpha',0.01);
[test_type_a(2,1) test_type_a(2,2)] = kstest2(burst_type_a(:,1),burst_type_a(:,2),'Alpha',0.05);
[test_type_a(3,1) test_type_a(3,2)] = kstest2(burst_type_a(:,1),burst_type_a(:,3),'Alpha',0.01);
[test_type_a(4,1) test_type_a(4,2)] = kstest2(burst_type_a(:,1),burst_type_a(:,3),'Alpha',0.05);
[test_type_a(5,1) test_type_a(5,2)] = kstest2(burst_type_a(:,2),burst_type_a(:,3),'Alpha',0.01);
[test_type_a(6,1) test_type_a(6,2)] = kstest2(burst_type_a(:,2),burst_type_a(:,3),'Alpha',0.05);

[test_type_b(1,1) test_type_b(1,2)] = kstest2(burst_type_b(:,1),burst_type_b(:,2),'Alpha',0.01);
[test_type_b(2,1) test_type_b(2,2)] = kstest2(burst_type_b(:,1),burst_type_b(:,2),'Alpha',0.05);
[test_type_b(3,1) test_type_b(3,2)] = kstest2(burst_type_b(:,1),burst_type_b(:,3),'Alpha',0.01);
[test_type_b(4,1) test_type_b(4,2)] = kstest2(burst_type_b(:,1),burst_type_b(:,3),'Alpha',0.05);
[test_type_b(5,1) test_type_b(5,2)] = kstest2(burst_type_b(:,2),burst_type_b(:,3),'Alpha',0.01);
[test_type_b(6,1) test_type_b(6,2)] = kstest2(burst_type_b(:,2),burst_type_b(:,3),'Alpha',0.05);


burst_type = [burst_type_z, burst_type_a, burst_type_b];


figure0 = figure;
set(gcf,'color','w')
subplot(1,3,1)
T = bplot(burst_type_z);
set(gca,'XTick',1:3,'XTickLabel',{'Before','Crisis','After'})

for i=1:3
    if test_type_z(2*(i-1) + 1 + alpha_ind,1) > 0
        sig_z(i,1) = test_type_z(2*(i-1) + 1 + alpha_ind,2);
    else
        sig_z(i,1) = nan;
    end
    if test_type_a(2*(i-1) + 1 + alpha_ind,1) > 0
        sig_a(i,1) = test_type_a(2*(i-1) + 1 + alpha_ind,2);
    else
        sig_a(i,1) = nan;
    end
    if test_type_b(2*(i-1) + 1 + alpha_ind,1) > 0
        sig_b(i,1) = test_type_b(2*(i-1) + 1 + alpha_ind,2);
    else
        sig_b(i,1) = nan;
    end
end

H = sigstar({[1,2], [1,3], [2,3]},[sig_z(1,1), sig_z(2,1), sig_z(3,1)]);
set(H,'color','k')
% ylim([.25 .55])
ylim([.05 .65])
xlabel('Period')
ylabel('Burstiness')
% ylabel('Memory Coefficient')
title('Type Zero')
% legend(T(1:4))
% legend boxoff

subplot(1,3,2)
T = bplot(burst_type_a);
set(gca,'XTick',1:3,'XTickLabel',{'Before','Crisis','After'})
% H = sigstar({[1,2], [1,3], [2,3]},[nan, test_type_a(3,2), nan]);
H = sigstar({[1,2], [1,3], [2,3]},[sig_a(1,1), sig_a(2,1), sig_a(3,1)]);
set(H,'color','k')
% ylim([.2 .48])
ylim([.05 .65])
xlabel('Period')
ylabel('Burstiness')
% ylabel('Memory Coefficient')
title('Type A')
% legend(T(1:4))
% legend boxoff

subplot(1,3,3)
T = bplot(burst_type_b);
set(gca,'XTick',1:3,'XTickLabel',{'Before','Crisis','After'})
H = sigstar({[1,2], [1,3], [2,3]},[sig_b(1,1), sig_b(2,1), sig_b(3,1)]);
set(H,'color','k')
% ylim([.07 .52])
ylim([.05 .65])
xlabel('Period')
ylabel('Burstiness')
% ylabel('Memory Coefficient')
title('Type B')
% legend(T(1:4))
% legend boxoff


dx=0.06;
dy=0.06;
x=(1-4*dx)/2.9;
y=(1-4*dy)/1;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure0,'Type','axes');
set(AxesHandle(3),'Position',[2*dx,3*dy,x,y]);
set(AxesHandle(2),'Position',[x+2.6*dx,3*dy,x,y]);
set(AxesHandle(1),'Position',[2*x+3.5*dx,3*dy,x,y]);
set(gcf,'Position',[100, 100, 1000, 500])


%% statistics for 30 days x 30 firms

addpath('/home/minyoung/data_minyoung/Research/Matlab_function/boxplot')
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/sigstar')

burst_b = [];
burst_c = [];
burst_a = [];
for i=1:30
    burst_b = [burst_b; mem_coef{i,1}];
    burst_c = [burst_c; mem_coef{i+30,1}];
    burst_a = [burst_a; mem_coef{i+130,1}];
end


burst_type_z(:,1) = burst_b(:,1);
burst_type_z(:,2) = burst_c(:,1);
burst_type_z(:,3) = burst_a(:,1);

burst_type_a(:,1) = burst_b(:,2);
burst_type_a(:,2) = burst_c(:,2);
burst_type_a(:,3) = burst_a(:,2);

burst_type_b(:,1) = burst_b(:,3);
burst_type_b(:,2) = burst_c(:,3);
burst_type_b(:,3) = burst_a(:,3);




test_b(1,1) = kstest2(burst_b(:,1),burst_b(:,2),'Alpha',0.01);
test_b(2,1) = kstest2(burst_b(:,1),burst_b(:,2),'Alpha',0.05);
test_b(3,1) = kstest2(burst_b(:,1),burst_b(:,3),'Alpha',0.01);
test_b(4,1) = kstest2(burst_b(:,1),burst_b(:,3),'Alpha',0.05);
test_b(5,1) = kstest2(burst_b(:,2),burst_b(:,3),'Alpha',0.01);
test_b(6,1) = kstest2(burst_b(:,2),burst_b(:,3),'Alpha',0.05);

test_c(1,1) = kstest2(burst_c(:,1),burst_c(:,2),'Alpha',0.01);
test_c(2,1) = kstest2(burst_c(:,1),burst_c(:,2),'Alpha',0.05);
test_c(3,1) = kstest2(burst_c(:,1),burst_c(:,3),'Alpha',0.01);
test_c(4,1) = kstest2(burst_c(:,1),burst_c(:,3),'Alpha',0.05);
test_c(5,1) = kstest2(burst_c(:,2),burst_c(:,3),'Alpha',0.01);
test_c(6,1) = kstest2(burst_c(:,2),burst_c(:,3),'Alpha',0.05);

test_a(1,1) = kstest2(burst_a(:,1),burst_a(:,2),'Alpha',0.01);
test_a(2,1) = kstest2(burst_a(:,1),burst_a(:,2),'Alpha',0.05);
test_a(3,1) = kstest2(burst_a(:,1),burst_a(:,3),'Alpha',0.01);
test_a(4,1) = kstest2(burst_a(:,1),burst_a(:,3),'Alpha',0.05);
test_a(5,1) = kstest2(burst_a(:,2),burst_a(:,3),'Alpha',0.01);
test_a(6,1) = kstest2(burst_a(:,2),burst_a(:,3),'Alpha',0.05);

% 
% figure
% subplot(1,3,1)
% bplot(burst_b)
% set(gca,'XTick',1:3,'XTickLabel',{'Zero','A','B'})
% H = sigstar({[1,2], [1,3], [2,3]},[0.01, 0.01, 0.01]);
% set(H,'color','k')
% ylim([.1 .55])
% 
% subplot(1,3,2)
% bplot(burst_c)
% set(gca,'XTick',1:3,'XTickLabel',{'Zero','A','B'})
% H = sigstar({[1,2], [1,3], [2,3]},[0.01, 0.01, nan]);
% set(H,'color','k')
% ylim([.1 .55])
% 
% subplot(1,3,3)
% bplot(burst_a)
% set(gca,'XTick',1:3,'XTickLabel',{'Zero','A','B'})
% H = sigstar({[1,2], [1,3], [2,3]},[0.01, 0.01, 0.01]);
% set(H,'color','k')
% ylim([.1 .55])

test_type_z = zeros(6,2);
test_type_a = zeros(6,2);
test_type_b = zeros(6,2);

[test_type_z(1,1) test_type_z(1,2)] = kstest2(burst_type_z(:,1),burst_type_z(:,2),'Alpha',0.01);
[test_type_z(2,1) test_type_z(2,2)] = kstest2(burst_type_z(:,1),burst_type_z(:,2),'Alpha',0.05);
[test_type_z(3,1) test_type_z(3,2)] = kstest2(burst_type_z(:,1),burst_type_z(:,3),'Alpha',0.01);
[test_type_z(4,1) test_type_z(4,2)] = kstest2(burst_type_z(:,1),burst_type_z(:,3),'Alpha',0.05);
[test_type_z(5,1) test_type_z(5,2)] = kstest2(burst_type_z(:,2),burst_type_z(:,3),'Alpha',0.01);
[test_type_z(6,1) test_type_z(6,2)] = kstest2(burst_type_z(:,2),burst_type_z(:,3),'Alpha',0.05);

[test_type_a(1,1) test_type_a(1,2)] = kstest2(burst_type_a(:,1),burst_type_a(:,2),'Alpha',0.01);
[test_type_a(2,1) test_type_a(2,2)] = kstest2(burst_type_a(:,1),burst_type_a(:,2),'Alpha',0.05);
[test_type_a(3,1) test_type_a(3,2)] = kstest2(burst_type_a(:,1),burst_type_a(:,3),'Alpha',0.01);
[test_type_a(4,1) test_type_a(4,2)] = kstest2(burst_type_a(:,1),burst_type_a(:,3),'Alpha',0.05);
[test_type_a(5,1) test_type_a(5,2)] = kstest2(burst_type_a(:,2),burst_type_a(:,3),'Alpha',0.01);
[test_type_a(6,1) test_type_a(6,2)] = kstest2(burst_type_a(:,2),burst_type_a(:,3),'Alpha',0.05);

[test_type_b(1,1) test_type_b(1,2)] = kstest2(burst_type_b(:,1),burst_type_b(:,2),'Alpha',0.01);
[test_type_b(2,1) test_type_b(2,2)] = kst0est2(burst_type_b(:,1),burst_type_b(:,2),'Alpha',0.05);
[test_type_b(3,1) test_type_b(3,2)] = kstest2(burst_type_b(:,1),burst_type_b(:,3),'Alpha',0.01);
[test_type_b(4,1) test_type_b(4,2)] = kstest2(burst_type_b(:,1),burst_type_b(:,3),'Alpha',0.05);
[test_type_b(5,1) test_type_b(5,2)] = kstest2(burst_type_b(:,2),burst_type_b(:,3),'Alpha',0.01);
[test_type_b(6,1) test_type_b(6,2)] = kstest2(burst_type_b(:,2),burst_type_b(:,3),'Alpha',0.05);


burst_type = [burst_type_z, burst_type_a, burst_type_b];


figure
set(gcf,'color','w')
subplot(1,3,1)
T = bplot(burst_type_z);
set(gca,'XTick',1:3,'XTickLabel',{'Before','Crisis','After'})
H = sigstar({[1,2], [1,3], [2,3]},[nan, nan, nan]);
set(H,'color','k')
% ylim([.25 .55])
% ylim([.05 .55])
xlabel('Period')
% ylabel('Burstiness')
ylabel('Memory Coefficient')
title('Type Zero')
legend(T(1:4))
legend boxoff

subplot(1,3,2)
T = bplot(burst_type_a);
set(gca,'XTick',1:3,'XTickLabel',{'Before','Crisis','After'})
% H = sigstar({[1,2], [1,3], [2,3]},[nan, test_type_a(3,2), nan]);
H = sigstar({[1,2], [1,3], [2,3]},[nan, nan, nan]);
set(H,'color','k')
% ylim([.2 .48])
% ylim([.05 .55])
xlabel('Period')
% ylabel('Burstiness')
ylabel('Memory Coefficient')
title('Type A')
legend(T(1:4))
legend boxoff

subplot(1,3,3)
T = bplot(burst_type_b);
set(gca,'XTick',1:3,'XTickLabel',{'Before','Crisis','After'})
H = sigstar({[1,2], [1,3], [2,3]},[test_type_b(1,2), nan, test_type_b(5,2)]);
set(H,'color','k')
% ylim([.07 .52])
xlabel('Period')
% ylabel('Burstiness')
ylabel('Memory Coefficient')
title('Type B')
legend(T(1:4))
legend boxoff


%%

addpath('/home/minyoung/data_minyoung/Research/Matlab_function/linspecer')


% burst & memory for plus
figure1 = figure(1);
set(gcf,'color','w')

subplot(1,2,1)
% plot(mean_burst_minus(:,1),'-v','color',[0 0 0])
% hold on
% plot(mean_burst_minus(:,2),'-o','color',[1 0 0])
% hold on
% plot(mean_burst_minus(:,3),'-o','color',[0 0 1])
hold on
% h = plot(mean_burst_plus(:,1),'-^','color',[0 0 0]);
% hold on
% plot(mean_burst_plus(:,2),'-o','markerfacecolor',[1 0 0],'color',[1 0 0])
% hold on
% plot(mean_burst_plus(:,3),'-o','color',[0 0 1])
% hold on

N = 3;
C = linspecer(N);
h1 = plot(mean_burst_plus(:,1),'-','color',C(3,:),'Linewidth',2);
hold on
h2 = plot(mean_burst_plus(:,2),'-','color',C(1,:),'Linewidth',2);
hold on
h3 = plot(mean_burst_plus(:,3),'-','color',C(2,:),'Linewidth',2);
hold on
plot([22 22],[.15 .46],'k--')
hold on
ylim([.15 .46])

date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
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
ax = ancestor(h3,'axes');
xrule = ax.XAxis;
xrule.FontSize = 8;
% ylim([0.14 0.39])

xtickangle(35)
xlim([1 160])
xlabel('Date')
ylabel('B_1')
lh = legend([h1 h2 h3],{'Type Zero, Positive','Type A, Positive','Type B, Positive'});
legend boxoff
set(lh, 'Position', [-.029 .26 .5 .01])
% legend('Zero Minus','A Minus','B Minus','Zero Plus','A Plus','B Plus')
% grid on
% grid minor
box on
text(-0.1,1.1,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')


subplot(1,2,2)
% plot(mean_mem_minus(:,1),'-v','color',[0 0 0])
% hold on
% plot(mean_mem_minus(:,2),'-o','color',[1 0 0])
% hold on
% plot(mean_mem_minus(:,3),'-o','color',[0 0 1])
hold on
h1 = plot(mean_mem_plus(:,1),'-','color',C(3,:),'Linewidth',2);
hold on
h2 = plot(mean_mem_plus(:,2),'-','color',C(1,:),'Linewidth',2);
hold on
h3 = plot(mean_mem_plus(:,3),'-','color',C(2,:),'Linewidth',2);
hold on
plot([22 22],[0 .25],'k--')
ylim([0 .25])
date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
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
ax = ancestor(h3,'axes');
xrule = ax.XAxis;
xrule.FontSize = 8;
xtickangle(35)
xlim([1 160])
xlabel('Date')
ylabel('M')
lh = legend([h1 h2 h3],{'Type Zero, Positive','Type A, Positive','Type B, Positive'});
legend boxoff
set(lh, 'Position', [.457 .26 .5 .01])
% legend('Zero Minus','A Minus','B Minus','Zero Plus','A Plus','B Plus')
% grid on
% grid minor
box on
text(-0.1,1.1,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')


dx=0.06;
dy=0.06;
x=(1-4*dx)/1.9;
y=(1-4*dy)/1.05;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure1,'Type','axes');
set(AxesHandle(2),'Position',[dx+0.01,3*dy,x,y]);
set(AxesHandle(1),'Position',[x+2.6*dx,3*dy,x,y]);
set(gcf,'Position',[100, 100, 900, 400])





% burst & memory for minus
figure2 = figure(2);
set(gcf,'color','w')

subplot(1,2,1)
% plot(mean_burst_minus(:,1),'-v','color',[0 0 0])
% hold on
% plot(mean_burst_minus(:,2),'-o','color',[1 0 0])
% hold on
% plot(mean_burst_minus(:,3),'-o','color',[0 0 1])
hold on
% h = plot(mean_burst_plus(:,1),'-^','color',[0 0 0]);
% hold on
% plot(mean_burst_plus(:,2),'-o','markerfacecolor',[1 0 0],'color',[1 0 0])
% hold on
% plot(mean_burst_plus(:,3),'-o','color',[0 0 1])
% hold on

N = 3;
C = linspecer(N);
h1 = plot(mean_burst_minus(:,1),'-','color',C(3,:),'Linewidth',2);
hold on
h2 = plot(mean_burst_minus(:,2),'-','color',C(1,:),'Linewidth',2);
hold on
h3 = plot(mean_burst_minus(:,3),'-','color',C(2,:),'Linewidth',2);
hold on
plot([22 22],[.12 .4],'k--')
hold on
ylim([.12 .4])

date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
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
ax = ancestor(h3,'axes');
xrule = ax.XAxis;
xrule.FontSize = 8;
% ylim([0.14 0.39])

xtickangle(35)
xlim([1 160])
xlabel('Date')
ylabel('B_1')
lh = legend([h1 h2 h3],{'Type Zero, Negative','Type A, Negative','Type B, Negative'});
legend boxoff
set(lh, 'Position', [-.029 .26 .5 .01])
% legend('Zero Minus','A Minus','B Minus','Zero Plus','A Plus','B Plus')
% grid on
% grid minor
box on
text(-0.1,1.1,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')


subplot(1,2,2)
% plot(mean_mem_minus(:,1),'-v','color',[0 0 0])
% hold on
% plot(mean_mem_minus(:,2),'-o','color',[1 0 0])
% hold on
% plot(mean_mem_minus(:,3),'-o','color',[0 0 1])
hold on
h1 = plot(mean_mem_minus(:,1),'-','color',C(3,:),'Linewidth',2);
hold on
h2 = plot(mean_mem_minus(:,2),'-','color',C(1,:),'Linewidth',2);
hold on
h3 = plot(mean_mem_minus(:,3),'-','color',C(2,:),'Linewidth',2);
hold on
plot([22 22],[-.1 .25],'k--')
ylim([-.01 .25])
date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
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
ax = ancestor(h3,'axes');
xrule = ax.XAxis;
xrule.FontSize = 8;
xtickangle(35)
xlim([1 160])
xlabel('Date')
ylabel('M')
lh = legend([h1 h2 h3],{'Type Zero, Negative','Type A, Negative','Type B, Negative'});
legend boxoff
set(lh, 'Position', [.457 .26 .5 .01])
% legend('Zero Minus','A Minus','B Minus','Zero Plus','A Plus','B Plus')
% grid on
% grid minor
box on
text(-0.1,1.1,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')


dx=0.06;
dy=0.06;
x=(1-4*dx)/1.9;
y=(1-4*dy)/1.05;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure2,'Type','axes');
set(AxesHandle(2),'Position',[dx+0.01,3*dy,x,y]);
set(AxesHandle(1),'Position',[x+2.6*dx,3*dy,x,y]);
set(gcf,'Position',[100, 100, 900, 400])



%%

addpath('/home/minyoung/data_minyoung/Research/Matlab_function/linspecer')


% % burst & memory for plus
% figure1 = figure(1);
% set(gcf,'color','w')
% 
% subplot(1,2,1)
% 
% N = 3;
% C = linspecer(N);
% h1 = plot(mean_burst_plus(:,1),'-','color',C(3,:),'Linewidth',1.5);
% hold on
% h2 = plot(mean_burst_plus(:,2),'-','color',C(1,:),'Linewidth',1.5);
% hold on
% h3 = plot(mean_burst_plus(:,3),'-','color',C(2,:),'Linewidth',1.5);
% hold on
% plot([22 22],[.115 .48],'k--')
% hold on
% ylim([.115 .48])
% 
% hold on
% ciplot(mean_burst_plus(:,1)-ci_burst_plus(:,1), mean_burst_plus(:,1)+ci_burst_plus(:,1), 1:160, C(3,:))
% hold on
% hold on
% ciplot(mean_burst_plus(:,2)-ci_burst_plus(:,2), mean_burst_plus(:,2)+ci_burst_plus(:,2), 1:160, C(1,:))
% hold on
% hold on
% ciplot(mean_burst_plus(:,3)-ci_burst_plus(:,3), mean_burst_plus(:,3)+ci_burst_plus(:,3), 1:160, C(2,:))
% hold on
% alpha(.5)
% 
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
% lh = legend([h1 h2 h3],{'Type Zero, Positive','Type A, Positive','Type B, Positive'});
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
% 
% N = 3;
% C = linspecer(N);
% h1 = plot(mean_mem_plus(:,1),'-','color',C(3,:),'Linewidth',1.5);
% hold on
% h2 = plot(mean_mem_plus(:,2),'-','color',C(1,:),'Linewidth',1.5);
% hold on
% h3 = plot(mean_mem_plus(:,3),'-','color',C(2,:),'Linewidth',1.5);
% hold on
% plot([22 22],[-.045 .26],'k--')
% hold on
% ylim([-.045 .26])
% 
% hold on
% ciplot(mean_mem_plus(:,1)-ci_mem_plus(:,1), mean_mem_plus(:,1)+ci_mem_plus(:,1), 1:160, C(3,:))
% hold on
% hold on
% ciplot(mean_mem_plus(:,2)-ci_mem_plus(:,2), mean_mem_plus(:,2)+ci_mem_plus(:,2), 1:160, C(1,:))
% hold on
% hold on
% ciplot(mean_mem_plus(:,3)-ci_mem_plus(:,3), mean_mem_plus(:,3)+ci_mem_plus(:,3), 1:160, C(2,:))
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
% lh = legend([h1 h2 h3],{'Type Zero, Positive','Type A, Positive','Type B, Positive'});
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
% AxesHandle=findobj(figure1,'Type','axes');
% set(AxesHandle(2),'Position',[dx+0.01,3*dy,x,y]);
% set(AxesHandle(1),'Position',[x+2.6*dx,3*dy,x,y]);
% set(gcf,'Position',[100, 100, 900, 400])





% burst & memory for minus
figure2 = figure(2);
set(gcf,'color','w')

subplot(1,2,1)

N = 3;
C = linspecer(N);
h1 = plot(mean_burst_minus(:,1),'-','color',C(3,:),'Linewidth',1.5);
hold on
h2 = plot(mean_burst_minus(:,2),'-','color',C(1,:),'Linewidth',1.5);
hold on
h3 = plot(mean_burst_minus(:,3),'-','color',C(2,:),'Linewidth',1.5);
hold on
h4 = plot(mean_burst_minus(:,4),'-','color',[1/3 1/3 1/3],'Linewidth',1.5);
hold on
plot([22 22],[.13 .51],'k--')
hold on
ylim([.13 .51])
% ylim([.13 .45])

hold on
ciplot(mean_burst_minus(:,1)-ci_burst_minus(:,1), mean_burst_minus(:,1)+ci_burst_minus(:,1), 1:160, C(3,:))
hold on
ciplot(mean_burst_minus(:,4)-ci_burst_minus(:,4), mean_burst_minus(:,4)+ci_burst_minus(:,4), 1:160, [5/8 5/8 5/8])
hold on
ciplot(mean_burst_minus(:,2)-ci_burst_minus(:,2), mean_burst_minus(:,2)+ci_burst_minus(:,2), 1:160, C(1,:))
hold on
ciplot(mean_burst_minus(:,3)-ci_burst_minus(:,3), mean_burst_minus(:,3)+ci_burst_minus(:,3), 1:160, C(2,:))
hold on
alpha(.5)

date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
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
ax = ancestor(h3,'axes');
xrule = ax.XAxis;
xrule.FontSize = 8;
% ylim([0.14 0.39])

xtickangle(35)
xlim([1 160])
xlabel('Date')
ylabel('B_1')
% lh = legend([h1 h2 h3],{'Type Zero','Type A','Type B'});
lh = legend([h1 h4 h2 h3],{'Type Zero','Type One', 'Type A','Type B'});
legend boxoff
set(lh, 'Position', [-.067 .26 .5 .01])
% legend('Zero Minus','A Minus','B Minus','Zero Plus','A Plus','B Plus')
% grid on
% grid minor
box on
text(-0.1,1.1,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')


subplot(1,2,2)
h1 = plot(mean_mem_minus(:,1),'-','color',C(3,:),'Linewidth',1.5);
hold on
h2 = plot(mean_mem_minus(:,2),'-','color',C(1,:),'Linewidth',1.5);
hold on
h3 = plot(mean_mem_minus(:,3),'-','color',C(2,:),'Linewidth',1.5);
hold on
h4 = plot(mean_mem_minus(:,4),'-','color',[1/3 1/3 1/3],'Linewidth',1.5);
hold on
plot([22 22],[0.02 .275],'k--')
hold on
ylim([.02 .275])

hold on
ciplot(mean_mem_minus(:,1)-ci_mem_minus(:,1), mean_mem_minus(:,1)+ci_mem_minus(:,1), 1:160, C(3,:))
hold on
ciplot(mean_mem_minus(:,4)-ci_mem_minus(:,4), mean_mem_minus(:,4)+ci_mem_minus(:,4), 1:160, [5/8 5/8 5/8])
hold on
ciplot(mean_mem_minus(:,2)-ci_mem_minus(:,2), mean_mem_minus(:,2)+ci_mem_minus(:,2), 1:160, C(1,:))
hold on
ciplot(mean_mem_minus(:,3)-ci_mem_minus(:,3), mean_mem_minus(:,3)+ci_mem_minus(:,3), 1:160, C(2,:))
hold on
alpha(.5)

date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
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
ax = ancestor(h3,'axes');
xrule = ax.XAxis;
xrule.FontSize = 8;
xtickangle(35)
xlim([1 160])
xlabel('Date')
ylabel('M')
% lh = legend([h1 h2 h3],{'Type Zero','Type A','Type B'});
lh = legend([h1 h4 h2 h3],{'Type Zero', 'Type One','Type A','Type B'});
legend boxoff
set(lh, 'Position', [.42 .26 .5 .01])
% legend('Zero Minus','A Minus','B Minus','Zero Plus','A Plus','B Plus')
% grid on
% grid minor
box on
text(-0.1,1.1,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')


dx=0.06;
dy=0.06;
x=(1-4*dx)/1.9;
y=(1-4*dy)/1.05;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure2,'Type','axes');
set(AxesHandle(2),'Position',[dx+0.01,3*dy,x,y]);
set(AxesHandle(1),'Position',[x+2.6*dx,3*dy,x,y]);
set(gcf,'Position',[100, 100, 900, 400])



%%

addpath('/home/minyoung/data_minyoung/Research/Matlab_function/linspecer')


% burst & memory for plus
figure1 = figure(1);
set(gcf,'color','w')

subplot(1,2,1)

N = 3;
C = linspecer(N);
h1 = plot(mean_burst_plus(:,1),'-','color',C(3,:),'Linewidth',1.5);
hold on
h2 = plot(mean_burst_plus(:,2),'-','color',C(1,:),'Linewidth',1.5);
hold on
h3 = plot(mean_burst_plus(:,3),'-','color',C(2,:),'Linewidth',1.5);
hold on
plot([22 22],[.115 .425],'k--')
hold on
ylim([-.01 .5])

hold on
ciplot(mean_burst_plus(:,1)-std_burst_plus(:,1), mean_burst_plus(:,1)+std_burst_plus(:,1), 1:160, C(3,:))
hold on
hold on
ciplot(mean_burst_plus(:,2)-std_burst_plus(:,2), mean_burst_plus(:,2)+std_burst_plus(:,2), 1:160, C(1,:))
hold on
hold on
ciplot(mean_burst_plus(:,3)-std_burst_plus(:,3), mean_burst_plus(:,3)+std_burst_plus(:,3), 1:160, C(2,:))
hold on
alpha(.5)


date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
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
ax = ancestor(h3,'axes');
xrule = ax.XAxis;
xrule.FontSize = 8;
% ylim([0.14 0.39])

xtickangle(35)
xlim([1 160])
xlabel('Date')
ylabel('B_1')
lh = legend([h1 h2 h3],{'Type Zero, Positive','Type A, Positive','Type B, Positive'});
legend boxoff
set(lh, 'Position', [-.029 .26 .5 .01])
% legend('Zero Minus','A Minus','B Minus','Zero Plus','A Plus','B Plus')
% grid on
% grid minor
box on
text(-0.1,1.1,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')


subplot(1,2,2)

N = 3;
C = linspecer(N);
h1 = plot(mean_mem_plus(:,1),'-','color',C(3,:),'Linewidth',1.5);
hold on
h2 = plot(mean_mem_plus(:,2),'-','color',C(1,:),'Linewidth',1.5);
hold on
h3 = plot(mean_mem_plus(:,3),'-','color',C(2,:),'Linewidth',1.5);
hold on
plot([22 22],[-.045 .25],'k--')
hold on
ylim([-.2 .35])

hold on
ciplot(mean_mem_plus(:,1)-std_mem_plus(:,1), mean_mem_plus(:,1)+std_mem_plus(:,1), 1:160, C(3,:))
hold on
hold on
ciplot(mean_mem_plus(:,2)-std_mem_plus(:,2), mean_mem_plus(:,2)+std_mem_plus(:,2), 1:160, C(1,:))
hold on
hold on
ciplot(mean_mem_plus(:,3)-std_mem_plus(:,3), mean_mem_plus(:,3)+std_mem_plus(:,3), 1:160, C(2,:))
hold on
alpha(.5)

date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
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
ax = ancestor(h3,'axes');
xrule = ax.XAxis;
xrule.FontSize = 8;
xtickangle(35)
xlim([1 160])
xlabel('Date')
ylabel('M')
lh = legend([h1 h2 h3],{'Type Zero, Positive','Type A, Positive','Type B, Positive'});
legend boxoff
set(lh, 'Position', [.457 .26 .5 .01])
% legend('Zero Minus','A Minus','B Minus','Zero Plus','A Plus','B Plus')
% grid on
% grid minor
box on
text(-0.1,1.1,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')


dx=0.06;
dy=0.06;
x=(1-4*dx)/1.9;
y=(1-4*dy)/1.05;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure1,'Type','axes');
set(AxesHandle(2),'Position',[dx+0.01,3*dy,x,y]);
set(AxesHandle(1),'Position',[x+2.6*dx,3*dy,x,y]);
set(gcf,'Position',[100, 100, 900, 400])





% burst & memory for minus
figure2 = figure(2);
set(gcf,'color','w')

subplot(1,2,1)

N = 3;
C = linspecer(N);
h1 = plot(mean_burst_minus(:,1),'-','color',C(3,:),'Linewidth',1.5);
hold on
h2 = plot(mean_burst_minus(:,2),'-','color',C(1,:),'Linewidth',1.5);
hold on
h3 = plot(mean_burst_minus(:,3),'-','color',C(2,:),'Linewidth',1.5);
hold on
plot([22 22],[.115 .425],'k--')
hold on
ylim([.0 .5])

hold on
ciplot(mean_burst_minus(:,1)-std_burst_minus(:,1), mean_burst_minus(:,1)+std_burst_minus(:,1), 1:160, C(3,:))
hold on
hold on
ciplot(mean_burst_minus(:,2)-std_burst_minus(:,2), mean_burst_minus(:,2)+std_burst_minus(:,2), 1:160, C(1,:))
hold on
hold on
ciplot(mean_burst_minus(:,3)-std_burst_minus(:,3), mean_burst_minus(:,3)+std_burst_minus(:,3), 1:160, C(2,:))
hold on
alpha(.5)

date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
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
ax = ancestor(h3,'axes');
xrule = ax.XAxis;
xrule.FontSize = 8;
% ylim([0.14 0.39])

xtickangle(35)
xlim([1 160])
xlabel('Date')
ylabel('B_1')
lh = legend([h1 h2 h3],{'Type Zero, Negative','Type A, Negative','Type B, Negative'});
legend boxoff
set(lh, 'Position', [-.029 .26 .5 .01])
% legend('Zero Minus','A Minus','B Minus','Zero Plus','A Plus','B Plus')
% grid on
% grid minor
box on
text(-0.1,1.1,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')


subplot(1,2,2)
h1 = plot(mean_mem_minus(:,1),'-','color',C(3,:),'Linewidth',1.5);
hold on
h2 = plot(mean_mem_minus(:,2),'-','color',C(1,:),'Linewidth',1.5);
hold on
h3 = plot(mean_mem_minus(:,3),'-','color',C(2,:),'Linewidth',1.5);
hold on
plot([22 22],[-.045 .25],'k--')
hold on
ylim([-.2 .34])

hold on
ciplot(mean_mem_minus(:,1)-std_mem_minus(:,1), mean_mem_minus(:,1)+std_mem_minus(:,1), 1:160, C(3,:))
hold on
hold on
ciplot(mean_mem_minus(:,2)-std_mem_minus(:,2), mean_mem_minus(:,2)+std_mem_minus(:,2), 1:160, C(1,:))
hold on
hold on
ciplot(mean_mem_minus(:,3)-std_mem_minus(:,3), mean_mem_minus(:,3)+std_mem_minus(:,3), 1:160, C(2,:))
hold on
alpha(.5)

date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
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
ax = ancestor(h3,'axes');
xrule = ax.XAxis;
xrule.FontSize = 8;
xtickangle(35)
xlim([1 160])
xlabel('Date')
ylabel('M')
lh = legend([h1 h2 h3],{'Type Zero, Negative','Type A, Negative','Type B, Negative'});
legend boxoff
set(lh, 'Position', [.457 .26 .5 .01])
% legend('Zero Minus','A Minus','B Minus','Zero Plus','A Plus','B Plus')
% grid on
% grid minor
box on
text(-0.1,1.1,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')


dx=0.06;
dy=0.06;
x=(1-4*dx)/1.9;
y=(1-4*dy)/1.05;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure2,'Type','axes');
set(AxesHandle(2),'Position',[dx+0.01,3*dy,x,y]);
set(AxesHandle(1),'Position',[x+2.6*dx,3*dy,x,y]);
set(gcf,'Position',[100, 100, 900, 400])

%%
figure
set(gcf,'color','w')

subplot(1,2,1)
% plot(mean_burst_minus(:,1),'-v','color',[0 0 0])
% hold on
% plot(mean_burst_minus(:,2),'-o','color',[1 0 0])
% hold on
% plot(mean_burst_minus(:,3),'-o','color',[0 0 1])
% hold on
hold on
plot(mean_burst_minus(:,1),'-^','color',[0 0 0])
hold on
plot(mean_burst_minus(:,2),'-o','markerfacecolor',[1 0 0],'color',[1 0 0])
hold on
plot(mean_burst_minus(:,3),'-o','markerfacecolor',[0 0 1],'color',[0 0 1])
hold on

eb = .5;
for i=1:size(mean_burst_minus,1)
    plot([i i],[mean_burst_minus(i,1)-std_burst_minus(i,1) mean_burst_minus(i,1)+std_burst_minus(i,1)],'k')
    hold on
    plot([i i],[mean_burst_minus(i,2)-std_burst_minus(i,2) mean_burst_minus(i,2)+std_burst_minus(i,2)],'r')
    hold on
    plot([i i],[mean_burst_minus(i,3)-std_burst_minus(i,3) mean_burst_minus(i,3)+std_burst_minus(i,3)],'b')
    hold on
    
    plot([i-eb i+eb],[mean_burst_minus(i,1)-std_burst_minus(i,1) mean_burst_minus(i,1)-std_burst_minus(i,1)],'k')
    hold on
    plot([i-eb i+eb],[mean_burst_minus(i,1)+std_burst_minus(i,1) mean_burst_minus(i,1)+std_burst_minus(i,1)],'k')
    hold on
    plot([i-eb i+eb],[mean_burst_minus(i,2)-std_burst_minus(i,2) mean_burst_minus(i,2)-std_burst_minus(i,2)],'r')
    hold on
    plot([i-eb i+eb],[mean_burst_minus(i,2)+std_burst_minus(i,2) mean_burst_minus(i,2)+std_burst_minus(i,2)],'r')
    hold on
    plot([i-eb i+eb],[mean_burst_minus(i,3)-std_burst_minus(i,3) mean_burst_minus(i,3)-std_burst_minus(i,3)],'b')
    hold on
    plot([i-eb i+eb],[mean_burst_minus(i,3)+std_burst_minus(i,3) mean_burst_minus(i,3)+std_burst_minus(i,3)],'b')
    hold on
end

% date_seq = load('/home/minyoung/data_minyoung/cm_code/test/date/date.txt');
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
xlabel('Date')
ylabel('Burstiness')
legend('Zero Minus','A Minus','B Minus')

subplot(1,2,2)
% plot(mean_mem_minus(:,1),'-v','color',[0 0 0])
% hold on
% plot(mean_mem_minus(:,2),'-o','color',[1 0 0])
% hold on
% plot(mean_mem_minus(:,3),'-o','color',[0 0 1])
% hold on
hold on
plot(mean_mem_minus(:,1),'-^','color',[0 0 0])
hold on
plot(mean_mem_minus(:,2),'-o','markerfacecolor',[1 0 0],'color',[1 0 0])
hold on
plot(mean_mem_minus(:,3),'-o','markerfacecolor',[0 0 1],'color',[0 0 1])


eb = .5;
for i=1:size(mean_mem_minus,1)
    plot([i i],[mean_mem_minus(i,1)-std_mem_minus(i,1) mean_mem_minus(i,1)+std_mem_minus(i,1)],'k')
    hold on
    plot([i i],[mean_mem_minus(i,2)-std_mem_minus(i,2) mean_mem_minus(i,2)+std_mem_minus(i,2)],'r')
    hold on
    plot([i i],[mean_mem_minus(i,3)-std_mem_minus(i,3) mean_mem_minus(i,3)+std_mem_minus(i,3)],'b')
    hold on
    plot([i-eb i+eb],[mean_mem_minus(i,1)-std_mem_minus(i,1) mean_mem_minus(i,1)-std_mem_minus(i,1)],'k')
    hold on
    plot([i-eb i+eb],[mean_mem_minus(i,1)+std_mem_minus(i,1) mean_mem_minus(i,1)+std_mem_minus(i,1)],'k')
    hold on
    plot([i-eb i+eb],[mean_mem_minus(i,2)-std_mem_minus(i,2) mean_mem_minus(i,2)-std_mem_minus(i,2)],'r')
    hold on
    plot([i-eb i+eb],[mean_mem_minus(i,2)+std_mem_minus(i,2) mean_mem_minus(i,2)+std_mem_minus(i,2)],'r')
    hold on
    plot([i-eb i+eb],[mean_mem_minus(i,3)-std_mem_minus(i,3) mean_mem_minus(i,3)-std_mem_minus(i,3)],'b')
    hold on
    plot([i-eb i+eb],[mean_mem_minus(i,3)+std_mem_minus(i,3) mean_mem_minus(i,3)+std_mem_minus(i,3)],'b')
    hold on
end
% set(gca,'xtick',date_first_ind)
% set(gca,'xticklabel',date_str_tick)
% xtickangle(45)
xlabel('Date')
ylabel('Memory Coefficient')
legend('Zero Minus','A Minus','B Minus')


%%
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

figure
set(gcf,'color','w')
subplot(3,1,1)
bar(ks_burst_minus(:,1))
title('KS test b/w Zero & A (Burstiness)')
subplot(3,1,2)
bar(ks_burst_minus(:,2))
title('KS test b/w Zero & B (Burstiness)')
subplot(3,1,3)
bar(ks_burst_minus(:,3))
title('KS test b/w A & B (Burstiness)')

%%

NI = 4;
win_size = 1;
mv_size = 1;
burst_ind = 1;

addpath('/home/minyoung/data_minyoung/Research/Matlab_function')
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/randraw')
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/boxplot')
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/linspecer')
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/granger_causality/')

load(sprintf('./burst_cal/burst_individual_1906_%d_%d_%d.mat',win_size,mv_size,burst_ind))
ABall = load(sprintf('./burst_cal/burst_result_ABall_%d_%d_%d.mat',win_size,mv_size,burst_ind));
mean_burst_minus_ABall = ABall.mean_burst_minus;
mean_burst_plus_ABall = ABall.mean_burst_plus;
mean_mem_minus_ABall = ABall.mean_mem_minus;
mean_mem_plus_ABall = ABall.mean_mem_plus;
clear ABall

date_seq = load('./date/date.txt');
market_stat_temp = load('/home/minyoung/data_minyoung/Research/Price_Impact/market_stat.mat');
market_stat_t = market_stat_temp.two_day; clear market_stat_temp
clear burst_stat_temp burst_stat
burst_stat_temp = load('/home/minyoung/data_minyoung/Research/Price_Impact/burst_exp_10000_5000.mat');
burst_stat(1,:) = burst_stat_temp.burst_mean;
burst_stat(2,:) = burst_stat_temp.burst_std;
burst_stat = burst_stat'; clear burst_stat_temp
FTSE_pre = load('ftse_daily.csv');
VIX_pre = load('vix.csv');
TEDspread_pre = load('ted_spread.csv');
load('snp500.mat');

market_stat(1:25,:) = market_stat_t(1:25,:);
market_stat(26,:) = [0, 0];
market_stat(27:169,:) = market_stat_t(26:end,:);
clear market_stat_t

VIX = zeros(length(date_seq),1);
TEDspread = zeros(length(date_seq),1);
SNP = zeros(length(date_seq),1);
for i=1:length(date_seq)
    for ii=1:length(VIX_pre(:,1))
        if date_seq(i,1) == VIX_pre(ii,1)
            VIX(i,1) = VIX_pre(ii,2);
        end
    end
    
    for ii=1:length(TEDspread_pre(:,1))
        if date_seq(i,1) == TEDspread_pre(ii,1)
            TEDspread(i,1) = TEDspread_pre(ii,2);
        end
    end
    
    for ii=1:length(snp(:,1))
        if date_seq(i,1) == snp(ii,1)
            SNP(i,1) = snp(ii,2);
        end
    end
end

% snp = zeros(size(snp2,1),2);
% for i=1:size(snp2,1)
%     snp(i,1) = str2num([snp2{i,1}{1}(1:4),snp2{i,1}{1}(6:7),snp2{i,1}{1}(9:10)]);
%     snp(i,2) = snp3(i,1);
% end

for i=1:length(VIX)
    if VIX(i,1) == 0
        VIX(i,1) = (VIX(i-1,1) + VIX(i+1,1))/2;
    end
    
    if TEDspread(i,1) == 0
        TEDspread(i,1) = (TEDspread(i-1,1) + TEDspread(i+1,1))/2;
    end
    
    if SNP(i,1) == 0
        SNP(i,1) = (SNP(i-1,1) + SNP(i+1,1))/2;
    end
end

FTSE = FTSE_pre(:,2);
clear VIX_pre TEDspread_pre i ii FTSE_pre snp

FTSE = moving_average(FTSE,win_size,mv_size,1);
VIX = moving_average(VIX,win_size,mv_size,1);
TEDspread = moving_average(TEDspread,win_size,mv_size,1);
market_stat = moving_sum(market_stat,win_size,mv_size,1);

sigma = 2;
end_point = 10000;
figure1 = figure(1);
set(gcf,'color','w')
% subplot(1,2,1)
% plot(burst_stat(2:end_point,1),'k-')
% hold on
% plot(burst_stat(2:end_point,1) - sigma * burst_stat(2:end_point,2),'b:')
% hold on
% ylabel('Burstiness')
% xlabel('Length')
% title('Poisson Process')
% plot(burst_stat(2:end_point,1) + sigma * burst_stat(2:end_point,2),'b:')
% xlim([0 3000])
% ylim([-.3 .3])
% subplot(1,2,2)
[p1 p2 p3] = plotyy(1:size(market_stat,1),market_stat(:,1),1:size(market_stat,1),market_stat(:,2));
xlabel('Days')
ylabel(p1(1),'Number of Market Order (Type A)')
ylabel(p1(2),'Number of Market Order (Type B)')
set(p1,{'ycolor'},{'k';'k'})
% title('Number of Market Order')
xlim(p1(1),[0 160])
xlim(p1(2),[0 160])
date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
    date_str{i,1} = num2str(date_seq(i,1));
end
date_first = [20080815;20080915;20081015;20081114;20081215;20090115;20090216;20090316];
for i=1:length(date_first)
    date_first_ind(i,1) = find(date_seq == date_first(i));
end
date_str = date_str';
date_str_tick = date_str(date_first_ind);

C = linspecer(3);
C(1,:) = C(3,:);

set(p2,'color',C(2,:),'linewidth',2)
set(p3,'color',C(3,:),'linewidth',2)
set(gca,'xtick',date_first_ind)
set(gca,'xticklabel',date_str_tick)
ax = ancestor(p2,'axes');
xrule = ax.XAxis;
xrule.FontSize = 8;
xtickangle(35)

dx=0.06;
dy=0.06;
x=(1-4*dx)/1;
y=(1-4*dy)/1;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure1,'Type','axes');
set(AxesHandle(1),'Position',[2*dx-.02,3*dy,x,y]);
set(gcf,'Position',[100, 100, 500, 400])
% date_str = cell(size(mean_burst_minus,1),1);
% for i=1:size(mean_burst_minus,1)
%     date_str{i,1} = num2str(date_seq(i,1));
% end
% date_first = [20080815;20080915;20081015;20081114;20081215;20090115;20090216;20090316];
% for i=1:length(date_first)
%     date_first_ind(i,1) = find(date_seq == date_first(i));
% end
% date_str = date_str';
% date_str_tick = date_str(date_first_ind)
% 
% set(gca,'xtick',date_first_ind)
% set(gca,'xticklabel',date_str_tick)
% xtickangle(45)


len = length(FTSE);

figure2 = figure(2);
set(gcf,'color','w')
subplot(1,2,1)

N = 3;
C = linspecer(N);
C1 = linspecer(10);
h = plot(mean_burst_plus(:,1),'-','color',C(3,:),'Linewidth',2);
hold on
plot(mean_burst_plus(:,2),'-','color',C(1,:),'Linewidth',2)
hold on
plot(mean_burst_plus(:,3),'-','color',C(2,:),'Linewidth',2)
hold on
% plot(mean_burst_plus_ABall(:,1),'-','color',C(4,:),'Linewidth',2)
% hold on
plot(mean_burst_plus_ABall(:,2),'-','color',C1(10,:),'Linewidth',2)
ylim([.09 .43])

date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
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
ax = ancestor(h,'axes');
xrule = ax.XAxis;
xrule.FontSize = 8;
% ylim([0.14 0.39])

xtickangle(35)
xlim([1 160])
xlabel('Date')
ylabel('B')
lh = legend('Zero Plus','A Plus','B Plus','All Plus');
legend boxoff
set(lh, 'Position', [-.07 .28 .5 .01])
% legend('Zero Minus','A Minus','B Minus','Zero Plus','A Plus','B Plus')
% grid on
% grid minor
box on
text(-0.1,1.1,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')

subplot(1,2,2)
% plot(mean_mem_minus(:,1),'-v','color',[0 0 0])
% hold on
% plot(mean_mem_minus(:,2),'-o','color',[1 0 0])
% hold on
% plot(mean_mem_minus(:,3),'-o','color',[0 0 1])
hold on
h = plot(mean_mem_plus(:,1),'-','color',C(3,:),'Linewidth',2);
hold on
plot(mean_mem_plus(:,2),'-','color',C(1,:),'Linewidth',2)
hold on
plot(mean_mem_plus(:,3),'-','color',C(2,:),'Linewidth',2)
hold on
plot(mean_mem_plus_ABall(:,2),'-','color',C1(10,:),'Linewidth',2)
date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
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
ax = ancestor(h,'axes');
xrule = ax.XAxis;
xrule.FontSize = 8;
xtickangle(35)
xlim([1 160])
xlabel('Date')
ylabel('M')
lh = legend('Zero Plus','A Plus','B Plus','All Plus');
legend boxoff
set(lh, 'Position', [.426 .28 .5 .01])
% legend('Zero Minus','A Minus','B Minus','Zero Plus','A Plus','B Plus')
% grid on
% grid minor
box on

dx=0.06;
dy=0.06;
x=(1-4*dx)/1.9;
y=(1-4*dy)/1;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure2,'Type','axes');
set(AxesHandle(2),'Position',[.01+dx,3*dy,x,y]);
set(AxesHandle(1),'Position',[x+2.6*dx,3*dy,x,y]);
set(gcf,'Position',[100, 100, 900, 400])


date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
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
xtickangle(45)

figure
subplot(2,4,3)
plot(mean_price)
% set(gca,'xtick',date_first_ind)
% set(gca,'xticklabel',date_str_tick)
% xtickangle(45)
xlabel('Date')
ylabel('Average Price')
xlim([0 len])

date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
    date_str{i,1} = num2str(date_seq(i,1));
end
date_first = [20080815;20080915;20081015;20081114;20081215;20090115;20090216;20090316];
for i=1:length(date_first)
    date_first_ind(i,1) = find(date_seq == date_first(i));
end
date_str = date_str';
date_str_tick = date_str(date_first_ind)

set(gca,'xtick',date_first_ind)
set(gca,'xticklabel',date_str_tick)
xtickangle(45)

subplot(2,4,4)
plot(mean_volatility)
% set(gca,'xtick',date_first_ind)
% set(gca,'xticklabel',date_str_tick)
% xtickangle(45)
xlabel('Date')
ylabel('Volatility')
xlim([0 len])

date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
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
xtickangle(45)

subplot(2,4,5)
plot(mean_trading_volume)
% set(gca,'xtick',date_first_ind)
% set(gca,'xticklabel',date_str_tick)
% xtickangle(45)
xlabel('Date')
ylabel('Trading Volume')
xlim([0 len])

date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
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
xtickangle(45)

subplot(2,4,6)
plot(VIX)
% set(gca,'xtick',date_first_ind)
% set(gca,'xticklabel',date_str_tick)
% xtickangle(45)
xlabel('Date')
ylabel('VIX')
xlim([0 len])

date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
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
xtickangle(45)

subplot(2,4,7)
plot(TEDspread)
% set(gca,'xtick',date_first_ind)
% set(gca,'xticklabel',date_str_tick)
% xtickangle(45)
xlabel('Date')
ylabel('TEDspread')
xlim([0 len])

date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
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
xtickangle(45)

subplot(2,4,8)
[p1 p2 p3] = plotyy(1:size(market_stat,1),market_stat(:,1),1:size(market_stat,1),market_stat(:,2));
xlabel('Time')
ylabel(p1(1),'# of Market Order, Type A')
ylabel(p1(2),'# of Market Order, Type B')
xlim(p1(1),[0 len])
xlim(p1(2),[0 len])

date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
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
xtickangle(45)
% 
%

series = zeros(size(VIX,1),27);
series(:,1:3) = mean_burst_minus(:,1:3);
series(:,4:6) = mean_burst_plus(:,1:3);
series(:,7:9) = mean_mem_minus(:,1:3);
series(:,10:12) = mean_mem_plus(:,1:3);
series(:,13) = mean_price;
series(:,14) = mean_trading_volume;
series(:,15) = mean_volatility;
series(:,16) = TEDspread;
series(:,17) = VIX;
series(:,18:19) = mean_burst_minus_ABall;
series(:,20:21) = mean_burst_plus_ABall;
series(:,22:23) = mean_mem_minus_ABall;
series(:,24:25) = mean_mem_plus_ABall;
% series(:,26) = market_stat(:,1);
% series(:,27) = market_stat(:,2);
series(:,26) = mean_burst_minus(:,4);
series(:,27) = mean_mem_minus(:,4);

if size(series,1) == 169
    [corr_series, corr_pval] = corr(series([1:25,27:169],:));
else
    [corr_series, corr_pval] = corr(series(:,:));
end

t_value_burstiness_minus = zeros(12,1);
ind = 0;
for i=1:3
    for ii=1:4
        if ii==4
            ind_ii = 19;
        else
            ind_ii = ii;
        end
        ind = ind + 1;
        tb = table(series(:,13+i),series(:,ind_ii),'VariableNames',{'var1','var2'});
        temp = fitlm(tb,'var1 ~ var2');
        t_value_burstiness_minus(ind,1) = table2array(temp.Coefficients(2,3));
        clear tb temp
    end
end
t_value_burstiness_plus = zeros(12,1);
ind = 0;
for i=1:3
    for ii=1:4
        if ii==4
            ind_ii = 21;
        else
            ind_ii = ii+3;
        end
        ind = ind + 1;
        tb = table(series(:,13+i),series(:,ind_ii),'VariableNames',{'var1','var2'});
        temp = fitlm(tb,'var1 ~ var2');
        t_value_burstiness_plus(ind,1) = table2array(temp.Coefficients(2,3));
        clear tb temp
    end
end
t_value_memory_minus = zeros(12,1);
ind = 0;
for i=1:3
    for ii=1:4
        if ii==4
            ind_ii = 23;
        else
            ind_ii = ii+6;
        end
        ind = ind + 1;
        tb = table(series(:,13+i),series(:,ind_ii),'VariableNames',{'var1','var2'});
        temp = fitlm(tb,'var1 ~ var2');
        t_value_memory_minus(ind,1) = table2array(temp.Coefficients(2,3));
        clear tb temp
    end
end
t_value_memory_plus = zeros(12,1);
ind = 0;
for i=1:3
    for ii=1:4
        if ii==4
            ind_ii = 25;
        else
            ind_ii = ii+6;
        end
        ind = ind + 1;
        tb = table(series(:,13+i),series(:,ind_ii),'VariableNames',{'var1','var2'});
        temp = fitlm(tb,'var1 ~ var2');
        t_value_memory_plus(ind,1) = table2array(temp.Coefficients(2,3));
        clear tb temp
    end
end
mean_t_value_burstiness = (t_value_burstiness_minus + t_value_burstiness_plus)/2;
mean_t_value_memory = (t_value_memory_minus + t_value_memory_plus)/2;


figure
set(gcf,'color','w')
subplot(3,4,1)
scatter(series(:,1),series(:,14),'.k')
subplot(3,4,2)
scatter(series(:,1),series(:,15),'.k')
subplot(3,4,3)
scatter(series(:,1),series(:,16),'.k')
subplot(3,4,4)
scatter(series(:,1),series(:,17),'.k')
subplot(3,4,5)
scatter(series(:,2),series(:,14),'.k')
subplot(3,4,6)
scatter(series(:,2),series(:,15),'.k')
subplot(3,4,7)
scatter(series(:,2),series(:,16),'.k')
subplot(3,4,8)
scatter(series(:,2),series(:,17),'.k')
subplot(3,4,9)
scatter(series(:,3),series(:,14),'.k')
subplot(3,4,10)
scatter(series(:,3),series(:,15),'.k')
subplot(3,4,11)
scatter(series(:,3),series(:,16),'.k')
subplot(3,4,12)
scatter(series(:,3),series(:,17),'.k')

% a1 = series(:,3);
% a2 = series(:,14);
% figure
% set(gcf,'color','w')
% subplot(1,2,1)
% dif_a1 = a1(2:end) - a1(1:end-1);
% dif_a2 = a2(2:end,:) - a2(1:end-1,:);
% for i=1:length(dif_a1)
%     ah = annotation('arrow','headstyle','cback3','headlength',5,'headwidth',5);
%     set(ah,'parent',gca);
%     set(ah,'position',[10*a1(i,1),10*a2(i,1),10*dif_a1(i,1),10*dif_a2(i,1)])
%     hold on
% %     quiver(a1(1:end-1,1),a2(1:end-1,3),dif_a1(:,1),dif_a2(:,1),0)
% end
% scatter(10*a1(1,1),10*a2(1,1),30,'o','filled','markerfacecolor',[0 0 1])
% hold on
% scatter(10*a1(end,1),10*a2(end,1),30,'o','markeredgecolor',[0 0 1])
% hold on
% scatter(10*a1(32,1),10*a2(32,1),30,'*','markeredgecolor',[0 0 1])
% hold on
% for i=1:15
%     scatter(10*a1(i*10,1),10*a2(i*10,1),30,'o','markeredgecolor',[0 0 1])
%     hold on
% end
% xlim([10*min(a1) 10*max(a1)])
% ylim([10*min(a2(:,1)) 10*max(a2(:,1))])
% % xlabel('Volatility')
% % ylabel('Burstiness')
% % title('Type Zero')
% % title(sprintf('Type B, %.2f, %.2f',corr(a1,a2(:,1)),corr(a1,mean_burst_plus(:,1))))
% % legend('Zero Minus','Zero Plus')
% subplot(1,2,2)
% for i=1:15
%     a(i,1) = 10*a1(i*10,1);
%     a(i,2) = 10*a2(i*10,1);
% end
% plot(a(:,1),a(:,2),'-*')



%%

market_stat(:,1) % length, A
market_stat(:,2) % length, B
burst_stat(:,1) % mean burst
burst_stat(:,2) % std burst
burst_len = 1:size(burst_stat,1);

len_series = size(mean_burst_minus,1);

for i=1:len_series
    len = find(burst_len == market_stat(i,1));
    for ii=1:3
        if isempty(len)
            len = size(burst_stat,1);
            burst_minus_nor(i,ii) = (mean_burst_minus(i,ii) - burst_stat(len,1)) / burst_stat(len,2);
            burst_plus_nor(i,ii) = (mean_burst_plus(i,ii) - burst_stat(len,1)) / burst_stat(len,2);
        else
            burst_minus_nor(i,ii) = (mean_burst_minus(i,ii) - burst_stat(len,1)) / burst_stat(len,2);
            burst_plus_nor(i,ii) = (mean_burst_plus(i,ii) - burst_stat(len,1)) / burst_stat(len,2);
        end
    end
end

figure
set(gcf,'color','w')
plot(burst_minus_nor(:,1),'-kv')
hold on
plot(burst_minus_nor(:,2),'-o','color',[1 0 0])
hold on
plot(burst_minus_nor(:,3),'-o','color',[0 0 1])
hold on
plot(burst_plus_nor(:,1),'k^')
hold on
plot(burst_plus_nor(:,2),'-o','color',[1 0 0],'markerfacecolor',[1 0 0])
hold on
plot(burst_plus_nor(:,3),'-o','color',[0 0 1],'markerfacecolor',[0 0 1])
xlabel('Days')
ylabel('Burstiness (std)')
legend('Zero Minus','A Minus','B Minus','Zero Plus','A Plus','B Plus')
xlim([0 length(date_seq)])



date_seq = load('./date/date.txt');
date_seq = date_seq(win_size:end);
date_str = cell(size(mean_burst_minus,1),1);
for i=1:size(mean_burst_minus,1)
    date_str{i,1} = num2str(date_seq(i,1));
end
date_first = [20080815;20080915;20081015;20081114;20081215;20090115;20090216;20090316];
for i=1:length(date_first)
    date_first_ind(i,1) = find(date_seq == date_first(i));
end
date_str = date_str';
date_str_tick = date_str(date_first_ind)

set(gca,'xtick',date_first_ind)
set(gca,'xticklabel',date_str_tick)
xtickangle(45)