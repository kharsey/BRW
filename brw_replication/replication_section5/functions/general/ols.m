function [b0,b1,omega,u] = ols(y,x,addconstant)
% function ols runs ols regression of y = alpha + beta * x + error
% INPUTS: y(T*Ny), x(T*Nx): data
%         addconstant: dummy for if we want to add a constant to x to
%               calculate intercept. 
%               = 0(default) if x already contains a column of ones
%               = 1 if x doesn't have a column of ones
% OUTPUT:   b0: intercept
%           b1: slope
%           omega: variance-covariance matrix


if nargin>2 && addconstant == 1
    T = length(y);
    x = [ones(T,1),x];
end
para = x\y;
b0 = para(1,:)';
b1 = para(2:end,:)';

u = y - x*para;
omega = var_OLS(u);