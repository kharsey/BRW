function var_u = var_OLS(u)
% var_OLS calculates the variance of the error terms for an OLS regression
% last updated 8/3/2010

[T,n] = size(u);
var_u = zeros(n);

for t = 1:T
    var_u = var_u + u(t,:)'*u(t,:);
end
var_u = var_u/T;