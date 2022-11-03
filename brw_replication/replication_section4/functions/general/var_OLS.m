function var_u = var_OLS(u)
[T,n] = size(u);
var_u = zeros(n);

for t = 1:T
    var_u = var_u + u(t,:)'*u(t,:);
end
var_u = var_u/T;