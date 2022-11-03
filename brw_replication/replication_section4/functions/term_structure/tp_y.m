function [yhat, yrn, ytp] = tp_y(mu, PHI, yields, mature, W, K1Q_X, rinfQ, sigma_cP, mat, dt)
% n: maturity of the forward rate of interest
Y_1 = yields*W';
T = length(yields);

mats = mature/dt;
iota = ones(3,1)*dt;
[PHIQ,muQ,delta1,delta0] = JSZ_transform(K1Q_X, rinfQ*dt, sigma_cP*dt,mats,W,iota); % should equal to K0Q_cP, K1Q_cP, rho0_cP, rho1_cP

[B,A] = ABdiff(mats,delta1,PHIQ,delta0,muQ,sigma_cP*dt);
[Brn,Arn] = ABdiff(mats,delta1,PHI,delta0,mu*dt,sigma_cP*dt);

Yhat = repmat(A',T,1)*12 + Y_1* B';
yhat = Yhat(:, mats==mat);

%Et_nminus1 = Et_Y1tplusn(mu, PHI,Y_1 ,n-1);
%E_ft = delta0+(Et_nminus1)*delta1;

Yrn = repmat(Arn',T,1)*12 + Y_1* Brn';
yrn = Yrn(:, mats==mat);

ytp = yhat - yrn;