function [Gamma_hat, alpha_hat, Omega_hat] = estVAR(X, p, demean, intercept)
% function [Gamma_hat, alpha_hat, Omega_hat] = regressVAR(X)
% inputs:
% X: T*N data matrix
% p: order of the VAR (default 1)
% demean: demean the data before estimation? (default false)
% intercept: include an intercept in the regression? (default true)
%
% output:
% Gamma_hat : N*N
% alpha_hat : N*1
% Omega_hat : N*N
%
% X(t+1) = alpha + Gamma*X(t) + eps(t+1), cov(eps(t+1)) = Omega
%
% Compute the maximum likelihood estimates of Gamma, alpha, Omega
%
% NOTE: The MLE estimates of Gamma, alpha do not depend on Omega.
% That is, the argmax_{Gamma,alpha} [L(X|Gamma,alpha,Omega)] = f(X)
% So this function provides MLE of Gamma, alpha for a fixed Omega.

% parameters/flags
if nargin<4; intercept=true; end;
if nargin<3; demean=false; end;
if nargin<2; p=1; end;

[T,k] = size(X);

if (demean)
    X = X - ones(T,1)*mean(X);
end

if (p==1)

    Yt = X(1:end-1,:);  % (T-1)*k
    Ytp1 = X(2:end,:);  % (T-1)*k

    Y = Ytp1.';  % k*(T-1) 
    if (intercept)
        Z = [ones(T-1,1), Yt].'; % (k+1)*(T-1)
    else
        Z = Yt.'; % k*(T-1)
    end

    A = Y*Z.'*inv(Z*Z.'); % k*(k+1) / k*k
    if (intercept)
        alpha_hat = A(:,1);
        Gamma_hat = A(:,2:end);
    else
        alpha_hat = NaN;
        Gamma_hat = A;
    end
    if nargout==3
        residuals = Ytp1 - (A*Z).'; % (T-1)*N
        Omega_hat = 1/(T-1)*(residuals')*residuals;
    end
elseif (p==2)
    Ytm2 = X(1:end-2,:);    % T-2 * k
    Ytm1 = X(2:end-1,:);    % T-2 * k
    Y = X(3:end,:)';      % k * T-2
    if (intercept)
        Z = [ones(T-2,1), Ytm1, Ytm2].'; % (2*k+1)*(T-2)
    else
        Z = [Ytm1, Ytm2].'; % 2k*(T-2)
    end
    A = Y*Z.'*inv(Z*Z.'); % k*(2k+1) / 2k*2k
    if (intercept)
        alpha_hat = A(:,1);
        Gamma_hat = A(:,2:end);
    else
        alpha_hat = NaN;
        Gamma_hat = A;
    end
    Gamma_hat = [Gamma_hat; [eye(k), zeros(k,k)]];
    if nargout==3
        residuals = Y - (A*Z).'; % (T-2)*N
        Omega_hat = 1/(T-2)*(residuals')*residuals;
    end
else
    error('not implemented for order>2');
end

