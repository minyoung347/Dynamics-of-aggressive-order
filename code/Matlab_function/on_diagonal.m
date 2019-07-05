function z = on_diagonal(x)

z = zeros(size(x,1),1);
for i = 1:size(x,1)
    z(i,1) = x(i,i);
end