function mcs = MCS_JSZ_pc_rest(parameters,zero_index,mature,W,mu_r,PHI_r,A_r,B_r,sig_r,H1,H1v,H2,W2)
% Calculate minimum chi-square statistic given structual, reduced form 
% parameters and Hessian matrices where there are zero restrictions on
% structural parameters.
% Notation: here r denotaes redueced form and s means structural parameters.


[lamQ,rinfQ,sig, sig_lam0,sig_lam1] = convert_v2ind_PoR(parameters,zero_index);
iota = ones(3,1)/12;
[PHIQ,muQ,delta1,delta0,penalty] = JSZ_transform(lamQ,rinfQ,sig,mature,W,iota);
if penalty > 100
    mcs = 1e50;
else
    [B_s,A_s] = ABdiff(mature,delta1,PHIQ,delta0,muQ,sig);
    mu_s = muQ + sig_lam0;
    PHI_s = PHIQ + sig_lam1;

    para1_diff = [mu_r';PHI_r'] - [mu_s';PHI_s'];
    F1 = para1_diff(:)'*H1*para1_diff(:);
    para1v_diff = sig_r*sig_r'- sig*sig';
    F1v = ltvec(para1v_diff)'*H1v*ltvec(para1v_diff);
    %F1v=0;
    para2_diff = [A_r';B_r']-[A_s';B_s']*W2';
    F2 = para2_diff(:)'*H2*para2_diff(:);
    mcs = F1+F1v+F2;
end