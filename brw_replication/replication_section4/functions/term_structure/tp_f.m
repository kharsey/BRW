function [f, frn, ftp] = tp_f(mu, PHI, yields, mature, W, K1Q_X, rinfQ, sigma_cP, n_from, n_to,dt)
% n: maturity of the forward rate of interest
Y_1 = yields*W';
T = length(yields);


mats = mature/dt;
iota = ones(3,1)*dt;
[PHIQ,muQ,delta1,delta0] = JSZ_transform(K1Q_X, rinfQ*dt, sigma_cP*dt,mats,W,iota); % should equal to K0Q_cP, K1Q_cP, rho0_cP, rho1_cP
%[B,A] = ABdiff(mats,delta1,PHIQ,delta0,muQ,sigma_cP/12); % SHOULD EQUAL TO BcP and AcP

mats = [n_from,n_to];
[B,A] = ABdiff(mats,delta1,PHIQ,delta0,muQ,sigma_cP*dt);
[Brn,Arn] = ABdiff(mats,delta1,PHI,delta0,mu*dt,sigma_cP*dt);

%{
% USING ONLY FUNCTION PROVIDED BY JSZ, USE WITH CAUTION.
mats = [n_from,n_to];
Sigma_cP = sigma_cP*sigma_cP';
K1Q_X = sort(K1Q_X,'descend');
K1Q_X = diag(K1Q_X-1);

[BcP, AcP, AX, BX] = jszLoadings(W, K1Q_X, rinfQ, Sigma_cP, mature, dt);
[K0Q_cP, K1Q_cP, rho0_cP, rho1_cP] = jszRotation(W, K1Q_X, rinfQ, dt, [], [], BX, AX);
[B,A] = gaussianDiscreteYieldLoadingsRecurrence(mats, K0Q_cP*dt, K1Q_cP,Sigma_cP*dt^2, rho0_cP*dt, rho1_cP);
[Brn,Arn] = gaussianDiscreteYieldLoadingsRecurrence(mats, mu*dt,PHI-eye(3), Sigma_cP*dt^2, rho0_cP*dt, rho1_cP);
%}

y_n = repmat(A',T,1)*12 + Y_1* B';
f = (y_n(:,2)*n_to - y_n(:,1)*n_from)/(n_to-n_from);

%Et_nminus1 = Et_Y1tplusn(mu, PHI,Y_1 ,n-1);
%E_ft = delta0+(Et_nminus1)*delta1;

yrn_n = repmat(Arn',T,1)*12 + Y_1* Brn';
frn = (yrn_n(:,2)*n_to - yrn_n(:,1)*n_from)/(n_to-n_from);

ftp = f - frn;