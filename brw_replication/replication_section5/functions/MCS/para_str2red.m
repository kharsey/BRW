function [A,B,penalty] = para_str2red(lamQ,rinfQ,sig,mature,W)
iota = ones(3,1)/12;
[PHIQ,muQ,delta1,delta0,penalty] = JSZ_transform(lamQ,rinfQ,sig,mature,W,iota);
if penalty > 100
    A = [];
    B = [];
else
    [B,A] = ABdiff(mature,delta1,PHIQ,delta0,muQ,sig);
end
