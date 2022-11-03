function [H, Hv] = information_matrix(x,Omega,addconstant)
% function information_matrix calculates the information matrix for OLS 
% estimates of an VAR
% INPUTS:   x, explanatory variables
%           Omega: variance-covariance matrix
%         addconstant: dummy for if we want to add a constant to x to
%               calculate intercept. 
%               = 0(default) if x already contains a column of ones
%               = 1(default) if x doesn't have a column of ones
% OUTPUTS:  H: information matrix for [b0, b1]
%           Hv: information matrix for Omega

if nargin>2 && addconstant == 1
    T = length(x);
    x = [ones(T,1),x];
end
T = length(x);
H = kron(inv(Omega),x'*x);
if nargout > 1
    d = duplication(size(Omega,1));
    Hv = T/2*d'*kron(inv(Omega),inv(Omega))*d;
end