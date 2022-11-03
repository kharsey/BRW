function mcs = MCS_JSZ_pc(parameters,mature,W,A_r,B_r,sig_r,H1v,H2,W2)
% Calculate minimum chi-square statistic given structual, reduced form 
% parameters and Hessian matrices
% Notation: here r denotaes redueced form and s means structural parameters.

[lamQ,rinfQ,sig] = convert_v2ind(parameters);
[A_s,B_s,penalty] = para_str2red(lamQ,rinfQ,sig,mature,W);

if penalty > 10
    mcs = 1e50;
else 
    para1v_diff = sig_r*sig_r'- sig*sig';
    F1v = ltvec(para1v_diff)'*H1v*ltvec(para1v_diff);
    %F1v=0;
    para2_diff = [A_r';B_r']-[A_s';B_s']*W2';
    F2 = para2_diff(:)'*H2*para2_diff(:);
    mcs = F1v+F2;
end