function [beta, alpha, Omega, log_llk] = reducedRankFreeInterceptRegress(Y, X, r, Omega)
% function [beta, alpha, Omega, log_llk] = reducedRankFreeInterceptRegress(Y, X, r, Omega)
% 
% Y : T*N
% X : T*M
% r : scalar, rankd of beta
% Omega : N*N
% 
% beta  : M*N
% alpha : 1*N
% Omega : N*N
% log_llk : T*1 minus log of the likelihood (including 2pi's)
%
% Regress Y(t,:) = X(t,:)*beta + alpha + eps(t), cov(eps) = Omega
%
%   under the assumption that beta has rank r, and alpha is unconstrained.
%  Note Letting Yhat(t,:) = Y(t,:) - Ybar
%               Xhat(t,:) = X(t,:) - Xbat
% for any Omega, beta, we can concentrate out alpha by rewriting the
% likelihood as 
% -N/2*log(2*pi) - .5*log(det(Omega)) 
%  - .5*sum((Yhat(t,:) - Xhat(t,:)*beta)*Omega^-1*(Yhat(t,:) - Xhat(t,:)*beta)
%  - T/2*(alpha - Ybar - beta'*Xbar).'(alpha - Ybar - beta'*Xbar)
% (cross terms cancel)
%
% The MLE of alpha is then alpha = Ybar = beta'*Xbar
%
% If Omega is not provided, likelihood is maximized over beta and Omega
%


[T N] = size(Y);
M = size(X,2);

Y0 = Y; % T*N
X0 = X; % T*M
Ybar = mean(Y); % 1*N
Xbar = mean(X); % 1*M
Y = Y - ones(T,1)*Ybar;
X = X - ones(T,1)*Xbar;

    
%M = size(X,2);

if nargin==2
    r = N;
end

if r==N
    beta = (X'*X)\(X'*Y); % M*N
else
    if nargin<4
        M00 = 1/(T) * Y'*Y;
        M0k = 1/(T) * Y'*X;
        Mkk = 1/(T) * X'*X;
        [V D] = eig(M0k'*(M00\M0k),Mkk);
        [A B] = sort(diag(D),'descend');
        b = V(:,B(1:r));
        %  b = b/b(1:r,1:r);
        a = M0k*b/(b'*Mkk*b);
        beta = (a*b')';
    else
        P = chol(X'*X);
        L = chol(Omega);
        betaols = (X'*X)\(X'*Y);
        [U S V] = svd(P*betaols/L);
        S(r+1:end,r+1:end) = 0;
        beta = P\U*S*V'*L;
    end
end

if nargout>2 && nargin<4
    Omega = 1/T*(Y-X*beta).'*(Y-X*beta);
end

if nargout>3
    log_llk = +N/2*log(2*pi) + .5*log(det(Omega)) + (sum((Y-X*beta).'.*(Omega\(Y-X*beta).'),1)).'; % T*1
end

alpha = Ybar - Xbar*beta;