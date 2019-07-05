function z = off_diagonal(x)

y = x.';
z = reshape(y(eye(size(x,1))'~=1),size(x,1)-1,size(x,1)).';
z = z(:);