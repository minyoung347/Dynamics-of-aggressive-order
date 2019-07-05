function d = kl_div(x1,x2,dir)
% Kullback_Leibler divergence, KLD
% x1 is [n x 1] original vector
% x2 is [n x 1] estimated vector
% if dir == 0, averaged kl_div, no direction
% if dir == 1, kl_div, direction

if size(x1,2) > 1 && size(x1,2) > 1 && size(x2,1) == 1 && size(x2,1) == 1
    x1 = x1';
    x2 = x2';
end

fig=figure; set(fig,'visible','off');
hist1 = histogram(x1,'Normalization','probability');
hist1_edge = hist1.BinEdges;
hist1_width = hist1.BinWidth;
close(fig);
fig=figure; set(fig,'visible','off');
hist2 = histogram(x2,'Normalization','probability');
hist2_edge = hist2.BinEdges;
hist2_width = hist2.BinWidth;
close(fig);
clear hist1 hist2

h_width = (hist1_width + hist2_width)/2;
h_minmax = [min(hist1_edge(1),hist2_edge(1)), max(hist1_edge(end),hist2_edge(end))];
h_edge = h_minmax(1) : h_width : h_minmax(2);
h_edge(length(h_edge)+1) = h_edge(end) + (h_edge(2) - h_edge(1));

fig=figure; set(fig,'visible','off');
hist1 = histogram(x1,h_edge,'Normalization','probability');
h1 = hist1.Values;
close(fig);
fig=figure; set(fig,'visible','off');
hist2 = histogram(x2,h_edge,'Normalization','probability');
h2 = hist2.Values;
close(fig);
clear hist1 hist2

d = zeros(size(h1));

% goodIdx = h1>0 & h2>0; 
% 
% d1 = sum(h1(goodIdx) .* log(h1(goodIdx) ./ h2(goodIdx)));
% d2 = sum(h2(goodIdx) .* log(h2(goodIdx) ./ h1(goodIdx)));

th = exp(-10);

h1(find(h1==0)) = th;
h2(find(h2==0)) = th;

d1 = sum(h1 .* log(h1 ./ h2));
d2 = sum(h2 .* log(h2 ./ h1));

if dir == 0
    d = (d1 + d2)/2;
else
    d = d1;
end

% x1 = randn(1000,1);
% a = histogram(x1);