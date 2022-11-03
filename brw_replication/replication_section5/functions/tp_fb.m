function [f, frn,ftp] = tp_f(mu, PHI,yields,mature, W,phiQ, rQ, sig,n_from,n_to,iota)
% n: maturity of the forward rate of interest
Y_1 = yields*W';
T = length(yields);

[PHIQ,muQ,delta1,delta0] = JSZ_transform(phiQ,rQ,sig,mature,W,iota);
mats = [n_from,n_to];
[B,A] = ABdiff(mats,delta1,PHIQ,delta0,muQ,sig);

y_n = repmat(A',T,1) + Y_1* B';
f = (y_n(:,2)*n_to - y_n(:,1)*n_from)/(n_to-n_from);

%Et_nminus1 = Et_Y1tplusn(mu, PHI,Y_1 ,n-1);
%E_ft = delta0+(Et_nminus1)*delta1;

[Brn,Arn] = ABdiff(mats,delta1,PHI,delta0,mu,sig);
yrn_n = repmat(Arn',T,1) + Y_1* Brn';
frn = (yrn_n(:,2)*n_to - yrn_n(:,1)*n_from)/(n_to-n_from);

ftp = f - frn;