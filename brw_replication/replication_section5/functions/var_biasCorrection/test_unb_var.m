%% load Wright data
start_quarter = 19901;
end_quarter = 20091;
wright_mats = [1,2,3,4,6,4*(2:10)];  % maturities in quarters
n_yields = length(wright_mats);

% load monthly yield file
yields = csvread('../../data/wright_us.csv');

% quarterly observations
yields = yields( mod(yields(:,2),3)==0,:);
quarter = yields(:,1)*10+yields(:,2)/3;

% sample period
yields = yields( (quarter>=start_quarter)&(quarter<=end_quarter),:);
T = length(yields);

Y = yields(:,2+wright_mats);
Y = Y/400;

% extract principal components
k = 3;
[V,D] = eigs(cov(Y),k);
V(:,[1,3]) = -V(:, [1,3]);
X = Y*V;


%% estimate VAR(1)
dbstop if error
%[Phi, mu, Omega] = estVAR(X, 1, true, false);
[Phi_tilde, mu_tilde, Omega_tilde] = est_unb_var(X, 2);%, true, 2, 2, 5000, true, 100000);
disp('largest root of bias-corrected Phi');
disp(max(eig(Phi_tilde)));

%% JSZ data
% break;
load('../jsz/sample_zeros.mat')
[V,D] = eigs(cov(yields),k);
X = yields*V;
[Phi_tilde, mu_tilde, Omega_tilde] = est_unb_var(X, 1);
[Phi_tilde, mu_tilde, Omega_tilde] = est_unb_var(X, 2);
 



