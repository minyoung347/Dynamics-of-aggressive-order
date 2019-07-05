function z = upper_diagonal(x)

temp = x;
z_temp = [];
for i = 1:size(temp,1)
    for ii = 1:size(temp,2)
        if i < ii
            z_temp = [z_temp; temp(i,ii)];
        end
    end
end

z = z_temp(:);