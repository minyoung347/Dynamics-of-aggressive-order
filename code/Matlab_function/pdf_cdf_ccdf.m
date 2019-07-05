%% linear binning using hist
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/randraw/')

a = randraw('pareto', [1, 2], [1 1e5]);
ss = length(a);

[pdf_temp,x] = hist(a,500);
resolution = x(2)-x(1);
pdf = pdf_temp/(ss*resolution);

cdf = cumsum(pdf*resolution);
ccdf = 1 - cdf;


figure
subplot(1,3,1)
plot(x,pdf)
set(gca,'xscale','log')
set(gca,'yscale','log')
subplot(1,3,2)
plot(x,cdf)
set(gca,'xscale','log')
set(gca,'yscale','log')
subplot(1,3,3)
plot(x,ccdf)
set(gca,'xscale','log')
set(gca,'yscale','log')




%% log_binning using hist
clear;clc;

addpath('/home/minyoung/data_minyoung/Research/Matlab_function/')
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/randraw/')

time_series = randraw('pareto', [1, 2], [1 1e5]);
% time_series = burst_dist_minus_ag{1,1};
len_time_series = length(time_series);
num_bin = 500;

% hist, linear1
[pdf_temp_linear_matf,x_linear_matf] = hist(time_series,num_bin);
df_linear_matf = x_linear_matf(2)-x_linear_matf(1);
% hist, linear2
[pdf_temp_linear_uf,x_linear_uf_trash] = hist(time_series,linspace(min(time_series),max(time_series),num_bin));
x_linear_uf_temp = linspace(min(time_series),max(time_series),num_bin+1);
df_linear_uf = x_linear_uf_temp(2:end) - x_linear_uf_temp(1:end-1);
x_linear_uf = (x_linear_uf_temp(1:end-1) + x_linear_uf_temp(2:end))/2;
% hist, log2
time_series = [1;5;11;99;101;1000];
num_bin = 3;
[pdf_temp_log_uf,x_log_uf_trash] = hist(time_series,10.^(linspace(log10(min(time_series)),log10(max(time_series)),num_bin+1)));
x_log_uf_temp = linspace(log10(min(time_series)),log10(max(time_series)),num_bin+1);
df_log_uf = 10.^(x_log_uf_temp(2:end)) - 10.^(x_log_uf_temp(1:end-1));
x_log_uf = 10.^((x_log_uf_temp(1:end-1) + x_log_uf_temp(2:end))/2);
% df_log_uf = log(x_log_uf(2)) - log(x_log_uf(1));
% x_log_uf = exp(linspace(log(min(time_series))+df_log_uf/2,log(max(time_series))-df_log_uf/2,num_bin));


% pdf, linear1
pdf_linear_matf = pdf_temp_linear_matf/(len_time_series*df_linear_matf);
% pdf, linear2
% pdf_linear_uf = pdf_temp_linear_uf/(len_time_series*df_linear_uf);
pdf_linear_uf = pdf_temp_linear_uf/sum(pdf_temp_linear_uf.*df_linear_uf);
% pdf, log2
pdf_log_uf = pdf_temp_log_uf/sum(pdf_temp_log_uf.*df_log_uf);


% cdf & ccdf, linear1
cdf_linear_matf = cumsum(pdf_linear_matf*df_linear_matf);
ccdf_linear_matf = 1 - cdf_linear_matf;
% cdf & ccdf, linear2
cdf_linear_uf = cumsum(pdf_linear_uf*df_linear_uf(1));
ccdf_linear_uf = 1 - cdf_linear_uf;
% cdf & ccdf, log2
cdf_log_uf = cumsum(pdf_log_uf.*df_log_uf);
ccdf_log_uf = 1 - cdf_log_uf;
% cdf_log_uf = cumsum(pdf_log_uf*df_log_uf);
% ccdf_log_uf = 1 - cdf_log_uf;



figure
subplot(1,3,1)
plot(x_linear_matf,pdf_linear_matf)
set(gca,'xscale','log')
set(gca,'yscale','log')
subplot(1,3,2)
plot(x_linear_matf,cdf_linear_matf)
set(gca,'xscale','log')
set(gca,'yscale','log')
subplot(1,3,3)
plot(x_linear_matf,ccdf_linear_matf)
set(gca,'xscale','log')
set(gca,'yscale','log')


figure
subplot(1,3,1)
plot(x_linear_uf,pdf_linear_uf)
set(gca,'xscale','log')
set(gca,'yscale','log')
subplot(1,3,2)
plot(x_linear_uf,cdf_linear_uf)
set(gca,'xscale','log')
set(gca,'yscale','log')
subplot(1,3,3)
plot(x_linear_uf,ccdf_linear_uf)
set(gca,'xscale','log')
set(gca,'yscale','log')


figure
subplot(1,3,1)
plot(x_log_uf,pdf_log_uf)
set(gca,'xscale','log')
set(gca,'yscale','log')
subplot(1,3,2)
plot(x_log_uf,cdf_log_uf)
set(gca,'xscale','log')
set(gca,'yscale','log')
subplot(1,3,3)
plot(x_log_uf,ccdf_log_uf)
set(gca,'xscale','log')
set(gca,'yscale','log')


%% linear_binning using histcounts

clear;clc;

addpath('/home/minyoung/data_minyoung/Research/Matlab_function/')
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/randraw/')

time_series = randraw('pareto', [1, 2], [1 1e5]);
% time_series = burst_dist_minus_ag{1,1};
len_time_series = length(time_series);


% hist, linear1
[h_value_temp, h_edge] = histcounts(time_series,'Normalization','probability');
x_mid = zeros(1,length(h_value_temp)+1);
for i=1:length(h_value_temp)
    mid_shift = abs(h_edge(i) - h_edge(i+1))/2;
    if i==1
        x_mid(i) = h_edge(i) - mid_shift;    
    end
        x_mid(i+1) = h_edge(i) + mid_shift;
end
h_value = [0, h_value_temp];

% cdf & ccdf, linear1
cdf_h = cumsum(h_value);
ccdf_h = 1 - cdf_h;


% % hist, linear1
% h_object = histogram(time_series,'Normalization','probability');
% h_num_bin = h_object.NumBins;
% h_edge = h_object.BinEdges;
% h_value_temp = h_object.Values;
% x_mid = zeros(1,h_num_bin+1);
% for i=1:h_num_bin
%     mid_shift = abs(h_edge(i) - h_edge(i+1))/2;
%     if i==1
%         x_mid(i) = h_edge(i) - mid_shift;    
%     end
%         x_mid(i+1) = h_edge(i) + mid_shift;
% end
% h_value = [0, h_value_temp];
% 
% % cdf & ccdf, linear1
% cdf_h = cumsum(h_value);
% ccdf_h = 1 - cdf_h;


figure
subplot(1,3,1)
plot(x_mid,h_value)
set(gca,'xscale','log')
set(gca,'yscale','log')
subplot(1,3,2)
plot(x_mid,cdf_h)
set(gca,'xscale','log')
set(gca,'yscale','log')
subplot(1,3,3)
plot(x_mid,ccdf_h)
set(gca,'xscale','log')
set(gca,'yscale','log')


%% log_binning using histogram

clear;clc;

addpath('/home/minyoung/data_minyoung/Research/Matlab_function/')
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/randraw/')

time_series = randraw('pareto', [1, 2], [1 1e5]);
% time_series = randraw('norm', [10, 1], 1, 1e5);
% time_series = 1:1e3;
log_time_series = log(time_series);
num_bin = 1000;

h_edge_log = linspace(min(log_time_series),max(log_time_series),num_bin+1);
h_df_log = h_edge_log(2:end) - h_edge_log(1:end-1);
h_count = zeros(1,num_bin);
for i=1:num_bin
    if i==num_bin
        h_count(i) = length(find(log_time_series >= h_edge_log(i) & log_time_series <= h_edge_log(i+1)));
    else
        h_count(i) = length(find(log_time_series >= h_edge_log(i) & log_time_series < h_edge_log(i+1)));
    end
end
h_edge = exp(1).^(h_edge_log);
h_df = h_edge(2:end) - h_edge(1:end-1);
h_mid = (h_edge(1:end-1) + h_edge(2:end))/2;

pdf_length = h_count./h_df;
pdf_length = pdf_length ./ length(time_series);
% pdf_numpoint = h_count./h_count;

cdf_length = cumsum(pdf_length .* h_df);
% cdf_numpoint = cumsum(pdf_numpoint);

ccdf_length = 1 - cdf_length;
% ccdf_numpoint = 1 - cdf_numpoint;
% 
% figure
% subplot(2,2,1)
% plot(h_mid)
% subplot(2,2,2)
% plot(h_count)
% subplot(2,2,3)
% plot(h_mid,h_count)
% subplot(2,2,4)
% plot(h_mid,h_count)
% set(gca,'xscale','log')
% set(gca,'yscale','log')


% figure
% subplot(1,3,1)
% plot(h_mid,pdf_length)
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% subplot(1,3,2)
% plot(h_mid,cdf_length)
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% subplot(1,3,3)
% plot(h_mid,ccdf_length)
% set(gca,'xscale','log')
% set(gca,'yscale','log')

[x_ccdf, y_ccdf] = ccdf(time_series);
figure
set(gcf,'color','w')
subplot(1,2,1)
plot(h_mid,ccdf_length)
set(gca,'xscale','log')
set(gca,'yscale','log')
ylim([10^-5 1])
subplot(1,2,2)
plot(x_ccdf,y_ccdf)
set(gca,'xscale','log')
set(gca,'yscale','log')




%% log_binning using histogram, ccdf = 1

clear;clc;

addpath('/home/minyoung/data_minyoung/Research/Matlab_function/')
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/randraw/')

% time_series = randraw('pareto', [1, 2], [1 1e5]);
time_series = randraw('norm', [10, 1], 1, 1e5);
len = length(time_series);
eps = 1/len;
% time_series = 1:1e3;
log_time_series = log(time_series);
num_bin = 1000;

h_edge_log_temp = linspace(min(log_time_series),max(log_time_series),num_bin);
temp = diff(h_edge_log_temp);
h_edge_log = [h_edge_log_temp(1)-temp(1), h_edge_log_temp];
h_df_log = h_edge_log(2:end) - h_edge_log(1:end-1);
h_count = zeros(1,num_bin);
for i=1:num_bin
    if i==num_bin
        h_count(i) = length(find(log_time_series >= h_edge_log(i) & log_time_series <= h_edge_log(i+1)));
    else
        h_count(i) = length(find(log_time_series >= h_edge_log(i) & log_time_series < h_edge_log(i+1)));
    end
end
h_edge = exp(1).^(h_edge_log);
h_df = h_edge(2:end) - h_edge(1:end-1);
h_mid = (h_edge(1:end-1) + h_edge(2:end))/2;

pdf_length = h_count./h_df;
pdf_length = pdf_length ./ length(time_series);
% pdf_numpoint = h_count./h_count;

cdf_length = cumsum(pdf_length .* h_df);
% cdf_numpoint = cumsum(pdf_numpoint);

ccdf_length = 1 - cdf_length;
ccdf_length = round(ccdf_length*len)/len;
% ccdf_numpoint = 1 - cdf_numpoint;
% 
% figure
% subplot(2,2,1)
% plot(h_mid)
% subplot(2,2,2)
% plot(h_count)
% subplot(2,2,3)
% plot(h_mid,h_count)
% subplot(2,2,4)
% plot(h_mid,h_count)
% set(gca,'xscale','log')
% set(gca,'yscale','log')


figure
subplot(1,3,1)
scatter(h_mid,pdf_length,'.')
xlim([min(h_mid) max(h_mid)])
set(gca,'xscale','log')
set(gca,'yscale','log')
subplot(1,3,2)
plot(h_mid,cdf_length)
xlim([min(h_mid) max(h_mid)])
set(gca,'xscale','log')
set(gca,'yscale','log')
subplot(1,3,3)
plot(h_mid,ccdf_length)
xlim([min(h_mid) max(h_mid)])
set(gca,'xscale','log')
set(gca,'yscale','log')

[x_ccdf, y_ccdf] = ccdf(time_series);
figure
set(gcf,'color','w')
subplot(1,2,1)
plot(h_mid,ccdf_length)
set(gca,'xscale','log')
set(gca,'yscale','log')
% ylim([10^-5 1])
subplot(1,2,2)
plot(x_ccdf,y_ccdf)
set(gca,'xscale','log')
set(gca,'yscale','log')



