clear
clc
addpath(genpath('functions'))

load('sample_zeros.mat')
forw = 2*yields(:,end)-yields(:,end-2);

mature = [6,12,24,36,60,84,120]/12;
J = length(mature);
dt = 1/12;
iota = ones(3,1)*dt;

options  =  optimset('fminunc');
options  =  optimset(options , 'TolFun'      , 1e-20);
options  =  optimset(options , 'TolX'        , 1e-15);
options  =  optimset(options , 'Display'     , 'on');
options  =  optimset(options , 'Diagnostics' , 'off');
options  =  optimset(options , 'LargeScale'  , 'off');
options  =  optimset(options , 'MaxFunEvals' , 10000) ;

cov_mat = cov(yields);
[eig_vec,eig_val] = eig(cov_mat);
eig_vec = eig_vec(:,end:-1:1);
eig_vec(:,1) = -eig_vec(:,1);

W = eig_vec(:,1:3)';
pc1 = yields*W';
T = length(pc1);

%% OLS

% P-dynamics
y1 = pc1(2:end,:);
x1 = [ones(T-1,1),pc1(1:end-1,:)];
para1 = x1\y1;
K0P_cP_ols  = para1(1,:)';
K1P_cP_ols = para1(2:end,:)';
u1t_star = y1 - x1*para1;
var_1 = var_OLS(u1t_star);
sigma_cP0 = chol(var_1)';

% Q-dynamics
% starting values
rinfQ = mean(yields(:,end));
K1Q_X = [.99,.98,.9];
parameters = convert_ind2v(sort(K1Q_X,'descend'),rinfQ,sigma_cP0);
% optimization
disp('optimization OLS');
[parameters,FVAL,EXITFLAG,OUTPUT] =  fminunc('jsz_like', parameters, options, yields, W, mature, dt);
[parameters,FVAL,EXITFLAG,OUTPUT] =  fminsearch('jsz_like', parameters, options, yields, W, mature, dt);
[parameters,FVAL,EXITFLAG,OUTPUT] =  fminunc('jsz_like', parameters, options, yields, W, mature, dt);

[K1Q_X_ols,rinfQ_ols,sigma_cP_ols] = convert_v2ind(parameters);

%% OLS standard errors
H1 = information_matrix(x1,var_1);
std_P_para = reshape(sqrt(diag(inv(H1))),4,3)';

[Vhat, Vhatrobust] = Vhat_MLE('jsz_like',parameters,yields,W, mature, dt);
std_para = sqrt(diag(Vhatrobust));
[std_K1Q_X,std_rinfQ,std_sigma_cP] = convert_v2ind(std_para);
% NOTE: we report mu*100, rinfQ*100, sigma_cP*100

%% output OLS results
% OUTPUT OF THE OLS PARAMETERS
disp('Table 1: 1st-3rd Columns')

%Parameters
K0P_cP_ols'*100
K1P_cP_ols

%Standard Errors
std_P_para(:,1)'*100
std_P_para(:,2:4)

sort(abs(eig(K1P_cP_ols)),'descend')'

[rinfQ_ols*100; std_rinfQ*100]

[K1Q_X_ols'; std_K1Q_X']

sigma_cP_ols*100
std_sigma_cP*100

% RMSE
mats = mature/dt;
[PHIQ,muQ,delta1,delta0] = JSZ_transform(K1Q_X_ols,rinfQ_ols*dt,sigma_cP_ols*dt,mats,W,iota);
[B,A] = ABdiff(mats,delta1,PHIQ,delta0,muQ,sigma_cP_ols*dt);
Yhat = repmat(A',T,1) + pc1*B';
rmse_ols = sqrt(sum(sum((yields - Yhat).^2))/T/J)*100;
disp(['RMSE = ', num2str(rmse_ols),' pts']);

%% Mean bias correction
% P-dynamics

mean_flag = true;

% inverse bootstrap
[Phi_tilde, mu_tilde, V_tilde] = est_unb_var(pc1, 1, mean_flag, 5000, 1000, 50, 1, 50000);
K0P_cP_mean  = mu_tilde;
K1P_cP_mean = Phi_tilde;
sigma_cP0 = chol(V_tilde)';

% Q-dynamics
% starting values
parameters = convert_ind2v(K1Q_X_ols,rinfQ_ols,sigma_cP0);
% optimization
disp('optimization MU');
[parameters] =  fminunc('jsz_like', parameters, options,yields,W, mature, dt, K0P_cP_mean, K1P_cP_mean-eye(3));
[parameters] =  fminsearch('jsz_like', parameters ,options,yields,W, mature, dt, K0P_cP_mean, K1P_cP_mean-eye(3));

[K1Q_X_mean,rinfQ_mean,sigma_cP_mean] = convert_v2ind(parameters);

% standard errors
H1 = information_matrix(x1,V_tilde);
std_P_para = reshape(sqrt(diag(inv(H1))),4,3)';

[Vhat, Vhatrobust] = Vhat_MLE('jsz_like', parameters, yields,W, mature, dt,K0P_cP_mean, K1P_cP_mean-eye(3));
std_para = sqrt(diag(Vhatrobust));
[std_K1Q_X,std_rinfQ,std_sigma_cP] = convert_v2ind(std_para);
% NOTE: we report K1Q_X, rinfQ*100, sigma_cP*100

%% OUPUT OF Mean-BC PARAMETERS
disp('Table 1: 4th-6th Columns')

disp('Mean-BC - P-dynamics');
K0P_cP_mean'*100
K1P_cP_mean
sort(abs(eig(K1P_cP_mean)),'descend')'
disp('Mean-BC - Q-dynamics');
rinfQ_mean*100
K1Q_X_mean'
disp('Mean-BC - Sigma');
sigma_cP_mean*100

disp('standard errors');
std_P_para(:,1)'*100
std_P_para(:,2:4)
std_rinfQ*100
std_K1Q_X'
std_sigma_cP*100

% RMSE
mats = mature/dt;
[PHIQ,muQ,delta1,delta0] = JSZ_transform(K1Q_X_mean,rinfQ_mean*dt,sigma_cP_mean*dt,mats,W,iota);
[B,A] = ABdiff(mats,delta1,PHIQ,delta0,muQ,sigma_cP_mean*dt);
Yhat = repmat(A',T,1) + pc1*B';
rmse_mean = sqrt(sum(sum((yields - Yhat).^2))/T/J)*100;
disp(['RMSE = ', num2str(rmse_mean),' pts']);


%% calculate forward rates
% 47-to-48 month forward rates
[f_ols, frn_ols,ftp_ols] = tp_f(K0P_cP_ols, K1P_cP_ols, yields, mature, W,K1Q_X_ols, rinfQ_ols, sigma_cP_ols,60,120,dt);
[f_mean, frn_mean,ftp_mean] = tp_f(K0P_cP_mean, K1P_cP_mean, yields, mature, W,K1Q_X_mean, rinfQ_mean, sigma_cP_mean,60,120,dt);


%% save all variables/estimates necessary for Monte Carlo study
estimates = ([K0P_cP_ols'; K1P_cP_ols; sigma_cP_ols; [rinfQ_ols,0,0]; K1Q_X_ols'] + ...
    [K0P_cP_mean'; K1P_cP_mean; sigma_cP_mean; [rinfQ_mean,0,0]; K1Q_X_mean'])/2;
estimates = [K0P_cP_mean'; K1P_cP_mean; sigma_cP_mean; [rinfQ_mean,0,0]; K1Q_X_mean'];

%% results
disp('Table 2: Output')

% largest eigenvalues
fprintf('$max(eig(Phi))$  ');
fprintf(' & %6.4f ', [max(abs(eig(K1P_cP_ols))), max(abs(eig(K1P_cP_mean)))]);
fprintf('\\\\ \n');

% half-life
max_irf = 12*40; 
irf_ols = irf_var1(K1P_cP_ols, max_irf);
hl_ols = find(irf_ols>.5,1,'last');
irf_mean = irf_var1(K1P_cP_mean, max_irf);
hl_mean = find(irf_mean>.5,1,'last');
fprintf('half-life        ');
fprintf(' &   %u ', [hl_ols, hl_mean]);
fprintf('\\\\ \n');

% CIR
tmp = inv(eye(3) - K1P_cP_ols);
cir_ols = tmp(1,1);
tmp = inv(eye(3) - K1P_cP_mean);
cir_mean = tmp(1,1);
fprintf('CIR              ');
fprintf(' & %6.1f ', [cir_ols, cir_mean]);
fprintf('\\\\ \n');

% IRF at 5y
fprintf('IRF at 5y        ');
fprintf(' & %6.2f ', [irf_ols(60), irf_mean(60)]);
fprintf('\\\\ \n');

% volatilities
fprintf('\\sigma(f_t^{60,120}) ');
fprintf(' & %4.3f ', 100*[std(f_ols), std(f_mean)]);
fprintf('\\\\ \n');

fprintf('\\sigma(\\tilde{f}_t^{60,120}) ');
fprintf(' & %4.3f ', 100*[std(frn_ols), std(frn_mean)]);
fprintf('\\\\ \n');

fprintf('\\sigma(\\tilde{f}_t^{60,120}) ');
fprintf(' & %4.3f ', 100*[std(ftp_ols), std(ftp_mean)]);
fprintf('\\\\ \n');

%% calculate forward rates
% 47-to-48 month forward rates
[f_ols, frn_ols,ftp_ols] = tp_f(K0P_cP_ols, K1P_cP_ols, yields, mature, W,K1Q_X_ols, rinfQ_ols, sigma_cP_ols,60,120,dt);
[f_mean, frn_mean,ftp_mean] = tp_f(K0P_cP_mean, K1P_cP_mean, yields, mature, W,K1Q_X_mean, rinfQ_mean, sigma_cP_mean,60,120,dt);

%% plot

[year month day] = datevec(dates);
yearIndex = [1; year(2:length(year))~=year(1:length(year)-1)]==1;
yearIndex = yearIndex.*(mod(year,5)==0);

figure(1);

set(gcf, 'Units', 'pixels')
set(gcf, 'Position', [100, 100, 700, 500])

min_yrange = -1;
max_yrange = ceil(100*max([forw;frn_ols;frn_mean;ftp_ols;ftp_mean]));
min_xrange = min(dates);
max_xrange = datenum(2010,1,1);

% recession data
startRec = datenum({'31-Jul-1990','31-Mar-2001','31-Dec-2007'});
endRec = datenum({'31-Mar-1991','30-Nov-2001','30-Jun-2009'});
nRec = 3;
y=[min_yrange max_yrange max_yrange min_yrange];
colorRec = [.7,.7,.7];

% plot risk-neutral rates
subplot('Position', [.08, .55, .9, .4]);
plot(dates, 100*forw, 'k-', 'Linewidth',2);
hold on;
for i=1:nRec
    x=[startRec(i) startRec(i) endRec(i) endRec(i)];
    h_fill = fill(x,y,colorRec);
    set(h_fill,'EdgeColor',colorRec);
end
set(gca, 'Layer', 'top')
h1 = plot(dates, 100*forw, 'k-', 'Linewidth',2);
h2 = plot(dates, 100*frn_ols, 'c-', 'Linewidth',2);
h3 = plot(dates, 100*frn_mean, 'b--', 'Linewidth',2);
ylabel('Percent');
axis([min_xrange max_xrange min_yrange max_yrange]);
line([min_xrange max_xrange],[0 0],'Color','k');
title('Risk-neutral rates');
legend([h1, h2, h3], {'forward rate', 'OLS', 'BC'});
set(gca, 'XTick', [dates(yearIndex==1);max_xrange]);
set(gca, 'XTickLabel', [year(yearIndex==1); 2010]);

% plot term premium
subplot('Position', [.08, .05, .9, .4]);
plot(dates, 100*ftp_ols, 'k-', 'Linewidth',2);
hold on;
for i=1:nRec
    x=[startRec(i) startRec(i) endRec(i) endRec(i)];
    h_fill = fill(x,y,colorRec);
    set(h_fill,'EdgeColor',colorRec);
end
set(gca, 'Layer', 'top')
h1 = plot(dates, 100*forw, 'k-', 'Linewidth',2);
h2 = plot(dates, 100*ftp_ols, 'c-', 'Linewidth',2);
h3 = plot(dates, 100*ftp_mean, 'b--', 'Linewidth',2);
ylabel('Percent');
axis([min_xrange max_xrange min_yrange max_yrange]);
line([min_xrange max_xrange],[0 0],'Color','k');
title('Forward term premia');
set(gca, 'XTick', [dates(yearIndex==1);max_xrange]);
set(gca, 'XTickLabel', [year(yearIndex==1); 2010]);
set(gcf, 'Color', 'w');

%% save decomposition (plot data) 
% [year,month] = datevec(dates);
% A = [year, month,forw, f_ols, frn_ols, ftp_ols, f_mean, frn_mean, ftp_mean];
% csvwrite('figure1_plotdata.csv',A);


