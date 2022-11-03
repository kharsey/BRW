function [K0Q_cP, K1Q_cP, rho0_cP, rho1_cP] = jszRotation(W, K1Q_X, rinfQ, dt, Sigma_cP, mats, BX, AX)
%function [K0Q_cP, K1Q_cP, rho0_cP, rho1_cP] = jszRotation(W, K1Q_X, rinfQ, dt, Sigma_cP, mats, BX, AX)
%
% Either provide (Sigma_cP, mats) or (BX, AX)           
%
% Inputs:
%   W          : N*J,      vector of portfolio weights to fit without error.
%   K1Q_X      : N*N
%   rinfQ      : scalar,   the long run mean under Q of the annualized short rate
%   dt         : scalar,   length of period in years
%   Sigma_cP   : N*N  covariance of innovations
%   mats       : 1*J,      maturities in years
%   BX         : N*J  (BX, AX) are optional (saves time)
%   AX         : 1*J
%
% Returns:
%   K0Q_cP : N*1
%   K1Q_cP : N*N
%   rho0_cP : scalar
%   rho1_cP : N*1
%
%
% r(t) = rho0_cP + rho1_cP'*cPt
%      = rinfQ + 1'*Xt
%      = 1 period discount rate (annualized)
%
% Under Q:
%   X(t+1) - X(t)   =          K1Q_X*X(t)  + eps_X(t+1),   cov(eps_X(t+1)) = Sigma_X
%   cP(t+1) - cP(t) = K0Q_cP + K1Q_cP*X(t) + eps_cP(t+1),  cov(eps_cP(t+1)) = Sigma_cP
%
% Where Sigma_X is chosen to match Sigma_cP 
%
% cPt = W*yt  (cPt N*1, W is N*J)
%     = W*AX' + W*BX'*Xt
%     = WAXp + WBXp*Xt
%
% Delta(cP) = WBXp*Delta(Xt)
%           = WBXp*(K1Q_X*Xt + sqrt(Sigma_X)*eps(t+1))
%           = WBXp*(K1Q_X)*(WBXp\(cPt - WAXp)) + sqrt(Sigma_cP)*eps(t+1)
%           = WBXp*(K1Q_X)/WBXp*cPt - WBXp*(K1Q_X)/WBXp*WAXp] + sqrt(Sigma_cP)*eps(t+1)
%
% rt = rinfQ + 1'*Xt  [annualized 1-period rate]
%    = rinfQ + 1'*(WBXp\(cPt - WAXp))
%    = [rinfQ - 1'*(WBXp\WAXp)] + ((WBXp)'1)'*cPt



if ~isempty(Sigma_cP)
    % Adjust K1Q_X in case we have near repeated root/zero eigenvalue
    K1Q_X = jszAdjustK1QX(K1Q_X);
    [BcP, AcP, AX, BX] = jszLoadings(W, K1Q_X, rinfQ, Sigma_cP, mats, dt);
end

N = size(K1Q_X,1);
WBXp = W*BX';
WAXp = W*AX';

K1Q_cP = WBXp*K1Q_X/WBXp;
K0Q_cP = - K1Q_cP*WAXp;

rho0_cP = rinfQ - ones(1,N)*(WBXp\WAXp);
rho1_cP = (WBXp)'\ones(N,1);
