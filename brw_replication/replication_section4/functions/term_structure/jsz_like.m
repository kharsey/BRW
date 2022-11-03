function [LLF,llk,penalty] = jsz_like(parameters,yields, W, mats, dt, K0P_cP, K1P_cP)
if max(parameters(1:3))>1
    LLF = 100000;
else
    K1Q_X = diag(parameters(1:3)-1);
    rinfQ = parameters(4);
    sig_p = veclt(parameters(5:end))/100;
    Sigma_cP = sig_p*sig_p';
    if nargin < 6
        [llk, penalty] = jszLLK(yields, W, K1Q_X, rinfQ, Sigma_cP, mats, dt);
    else
    [llk, penalty] = jszLLK(yields, W, K1Q_X, rinfQ, Sigma_cP, mats, dt, K0P_cP, K1P_cP);
    end
    LLF = sum(llk);
end