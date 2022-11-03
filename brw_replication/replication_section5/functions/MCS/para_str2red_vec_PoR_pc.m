function para_reduced = para_str2red_vec_PoR_pc(para_str,mature,W,W2,zero_index)
if nargin > 4
    [lamQ,rinfQ,sig, sig_lam0,sig_lam1] = convert_v2ind_PoR(para_str,zero_index);
else
    [lamQ,rinfQ,sig, sig_lam0,sig_lam1] = convert_v2ind_PoR(para_str);
end
iota = ones(3,1)/12;
[PHIQ,muQ,delta1,delta0,penalty] = JSZ_transform(lamQ,rinfQ,sig,mature,W,iota);
if penalty > 100
    A = [];
    B = [];
else
    [B,A] = ABdiff(mature,delta1,PHIQ,delta0,muQ,sig);
end

mu = muQ+sig_lam0; 
PHI = PHIQ + sig_lam1;

para_reduced1 = [mu';PHI'];
para_reduced1v = sig*sig';
para_reduced2 = [A';B']*W2';

para_reduced = [para_reduced1(:);ltvec(para_reduced1v);para_reduced2(:)];
