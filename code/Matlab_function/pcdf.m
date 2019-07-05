function [h_mid,pdf_length,cdf_length,ccdf_length] = pcdf(time_series,num_bin)
% [h_mid,pdf_length,cdf_length,ccdf_length] = pcdf(time_series,num_bin)
% calculate pdf, cdf, ccdf of the given time_series using log_binning
% num_bin is the number of log_bins
% h_mid is the middle points of log_bins, arithmetic mean (a+b)/2;
% only valid for positive values, time_series(find(time_series>0))

% example time_series
% addpath('/home/minyoung/data_minyoung/Research/Matlab_function/')
% addpath('/home/minyoung/data_minyoung/Research/Matlab_function/randraw/')
% time_series = 1:1e3;
% time_series = randraw('norm', [10, 1], 1, 1e5);
% time_series = randraw('pareto', [1, 2], [1 1e5]);
% num_bin = 1000; % for example

% only vaild for potitive values, time_series(find(time_series>0))
time_series = time_series(find(time_series>0));

len = length(time_series);
log_time_series = log(time_series);

h_edge_log_temp = linspace(min(log_time_series),max(log_time_series),num_bin); % n-1 bins of log_time_series, edge(n)
temp = diff(h_edge_log_temp); % difference between bins in log scale
h_edge_log = [h_edge_log_temp(1)-temp(1), h_edge_log_temp]; % n bins of log_time_series, edge(n+1), to include zero probability, ccdf=1
h_df_log = h_edge_log(2:end) - h_edge_log(1:end-1); % size of bins
% count the number of samples in each bin
h_count = zeros(1,num_bin);
for i=1:num_bin
    if i==num_bin
        h_count(i) = length(find(log_time_series >= h_edge_log(i) & log_time_series <= h_edge_log(i+1)));
    else
        h_count(i) = length(find(log_time_series >= h_edge_log(i) & log_time_series < h_edge_log(i+1)));
    end
end

h_edge = exp(1).^(h_edge_log); % edge of original time_series
h_df = h_edge(2:end) - h_edge(1:end-1); % difference between bins in original scale
h_mid = (h_edge(1:end-1) + h_edge(2:end))/2; % middle point of bins, arithmetic mean (a+b)/2;

% calculate pdf(probability density function)
pdf_length = h_count./h_df; % ex) 9/9, 90/90, 900/900 for uniform distribution in log binning
pdf_length = pdf_length ./ length(time_series);
% pdf_numpoint = h_count./h_count; % can't use pdf_numpoint in pdf

% calculate cdf(cumulative density function)
cdf_length = cumsum(pdf_length .* h_df);
% cdf_numpoint = cumsum(pdf_numpoint); % can't use pdf_numpoint in pdf

% calculate ccdf(complementary cumulative density function)
ccdf_length = 1 - cdf_length;
ccdf_length = round(ccdf_length*len)/len; % handling floating point
% ccdf_numpoint = 1 - cdf_numpoint; % can't use pdf_numpoint in pdf


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

% plot pdf, cdf, ccdf
% figure
% subplot(1,3,1)
% scatter(h_mid,pdf_length,'.')
% xlim([min(h_mid) max(h_mid)])
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% subplot(1,3,2)
% plot(h_mid,cdf_length)
% xlim([min(h_mid) max(h_mid)])
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% subplot(1,3,3)
% plot(h_mid,ccdf_length)
% xlim([min(h_mid) max(h_mid)])
% set(gca,'xscale','log')
% set(gca,'yscale','log')

% compare two ccdf
% [x_ccdf, y_ccdf] = ccdf(time_series);
% figure
% set(gcf,'color','w')
% subplot(1,2,1)
% plot(h_mid,ccdf_length)
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% subplot(1,2,2)
% plot(x_ccdf,y_ccdf)
% set(gca,'xscale','log')
% set(gca,'yscale','log')



