function [xcdf, yccdf]= ccdf(time_series)

time_series=time_series(find(time_series>0));

% figure
% set(gcf,'color','w')
[ycdf,xcdf] = cdfcalc(time_series);
xccdf = xcdf;
yccdf = 1-ycdf(1:end-1);
% if dash==1
%     plot(xcdf,yccdf,'--','LineWidth',1,'Color',color);
% elseif dash==0
%     plot(xcdf,yccdf,'LineWidth',1,'Color',color);
% end
% hold on
% % xlim([min(time_series)*.9, max(time_series)*1.2])
% set(gca,'xscale','log')
% set(gca,'yscale','log')