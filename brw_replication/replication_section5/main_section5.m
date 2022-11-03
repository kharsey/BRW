%==========================================================================
% Estimate unrestricted and restricted DTSM
% restrictions on sig*lambda_1
% chosen according to t-stats
% restrictions specified in convert_xxxx_PoR_rest.m
% @BRW, 2011-Nov
%==========================================================================

clear
clc
addpath(genpath('functions')) 

% GSW data
load('gsw_data_monthly.mat');
smpl = dates>datenum(1985,1,1) & dates<datenum(2012,1,1);
dates = dates(smpl);
yield_select = [1,2,3,5,7,10];  %(1:10);
J = length(yield_select);
yields = Y(smpl,yield_select)/1200;
mature = 12*yield_select;
forw = 2*yields(:,yield_select==10)-yields(:,yield_select==5); % five-by-ten-year forward rate

%% principal components and risk factors
cov_mat = cov(yields);   
[eig_vec,eig_val] = eig(cov_mat);
eig_vec = eig_vec(:,end:-1:1);

W = eig_vec(:,1:3)';
W2 = eye(J);

Y1 = yields*W';
Y2 = yields*W2';

%% FIRST STAGE -- regressions
T = length(Y1);

% first equation
[mu_ols_ur,PHI_ols_ur,var_1] = ols(Y1(2:end,:),Y1(1:end-1,:),1);
Sigma = chol(var_1)';

% second equation
[WA,WB,DUMMY,u2] = ols(Y2(2:end,:),Y1(2:end,:),1);
SIGe = sqrt(sum(sum((u2).^2))/(T-1)/size(u2,2));        % follow JSZ, the measurement errors are identically distributed across different matirities
var_2 = SIGe^2*eye(size(u2,2)); 

%------------------information matrix---------------
[H1, H1v] = information_matrix(Y1(1:end-1,:),var_1,1);
H2 = information_matrix(Y1(2:end,:),var_2,1);
std_P_para = reshape(sqrt(diag(inv(H1))),4,3)';

%% SECOND STAGE -- UNRESTRICTED OLS
options  =  optimset('fminunc');
options  =  optimset(options , 'TolFun'      , 1e-20);
options  =  optimset(options , 'TolX'        , 1e-15);
options  =  optimset(options , 'Display'     , 'on');
options  =  optimset(options , 'Diagnostics' , 'off');
options  =  optimset(options , 'LargeScale'  , 'off');
options  =  optimset(options , 'MaxFunEvals' , 50000) ;

% assign initial values
lamQ_start = [.99;.95;.9];
rinfQ_start = mean(yields(:,end));
sig_start = Sigma;
parameters = convert_ind2v(lamQ_start,rinfQ_start,sig_start);

disp('optimization -- OLS UR');
[parameters,FVAL,EXITFLAG,OUTPUT] =  fminunc('MCS_JSZ_pc',parameters ,options,mature,W,WA,WB,Sigma,H1v,H2,W2);

[lamQ_ols_ur,rinfQ_ols_ur,sig_ols_ur] = convert_v2ind(parameters);

%----------------standard errors for sig_lam representation----------------
iota = ones(3,1)/12;                        % iota = ones(3,1) in JSZ, but we are working with different scale here. Through I don't see any difference
[PHIQ_ols_ur,muQ_ols_ur] = JSZ_transform(lamQ_ols_ur,rinfQ_ols_ur,sig_ols_ur,mature,W,iota);
sig_lam0_ols_ur = mu_ols_ur - muQ_ols_ur;
sig_lam1_ols_ur = PHI_ols_ur - PHIQ_ols_ur;
R_PoR = blkdiag(H1,H1v,H2);
para_PoR = convert_ind2v_PoR(lamQ_ols_ur,rinfQ_ols_ur,sig_ols_ur, sig_lam0_ols_ur,sig_lam1_ols_ur);
GAMMA_PoR = numgrad('para_str2red_vec_PoR_pc',para_PoR,mature,W,W2);
variance_PoR = inv(GAMMA_PoR'*R_PoR*GAMMA_PoR);
para_std_PoR = sqrt(diag(variance_PoR));
[lamQ_ols_ur_std,rinfQ_ols_ur_std,sig_ols_ur_std, sig_lam0_ols_ur_std,sig_lam1_ols_ur_std] = convert_v2ind_PoR(para_std_PoR);

disp('t-stats of sig*lam1');
t_stats = sig_lam1_ols_ur ./ sig_lam1_ols_ur_std

% show RMSE for this specification
[PHIQ,muQ,delta1,delta0] = JSZ_transform(lamQ_ols_ur,rinfQ_ols_ur,sig_ols_ur,mature,W,iota);
[B,A] = ABdiff(mature,delta1,PHIQ,delta0,muQ,sig_ols_ur);
Yhat = repmat(A',T,1) + Y1*B';
rmse = sqrt(sum(sum((yields - Yhat).^2))/T/J)*1200;
disp(['RMSE = ', num2str(rmse),' pts']);

%% SECOND STAGE -- RESTRICTED OLS

% restrictions on sig*lam1 -- automatically chosen based on t-stats
zero_restr = ones(3,3);
zero_restr(abs(t_stats)<1) = nan;

zero_index = isnan(convert_ind2v_PoR(lamQ_ols_ur,rinfQ_ols_ur,sig_ols_ur, sig_lam0_ols_ur,zero_restr));
parameters = convert_ind2v_PoR(lamQ_ols_ur,rinfQ_ols_ur,sig_ols_ur, sig_lam0_ols_ur,sig_lam1_ols_ur,zero_index);

disp('optimization -- OLS R');
[parameters,FVAL,EXITFLAG,OUTPUT] =  fminunc('MCS_JSZ_pc_rest', parameters, options, zero_index, mature,W,mu_ols_ur,PHI_ols_ur,WA,WB,Sigma,H1,H1v,H2,W2);
[parameters,FVAL,EXITFLAG,OUTPUT] =  fminsearch('MCS_JSZ_pc_rest', parameters, options, zero_index, mature,W,mu_ols_ur,PHI_ols_ur,WA,WB,Sigma,H1,H1v,H2,W2);

[lamQ_ols_r,rinfQ_ols_r,sig_ols_r, sig_lam0_ols_r,sig_lam1_ols_r] = convert_v2ind_PoR(parameters,zero_index);
[PHIQ_ols_r,muQ_ols_r] = JSZ_transform(lamQ_ols_r,rinfQ_ols_r,sig_ols_r,mature,W,iota);
mu_ols_r = sig_lam0_ols_r + muQ_ols_r;
PHI_ols_r = sig_lam1_ols_r + PHIQ_ols_r;

%----------------standard errors for sig_lam representation----------------
R_PoR = blkdiag(H1,H1v,H2);
GAMMA_PoR = numgrad('para_str2red_vec_PoR_pc',parameters,mature,W,W2,zero_index);
variance_PoR = inv(GAMMA_PoR'*R_PoR*GAMMA_PoR);
para_std_PoR = sqrt(diag(variance_PoR));
[lamQ_ols_r_std,rinfQ_ols_r_std,sig_ols_r_std, sig_lam0_ols_r_std,sig_lam1_ols_r_std] = convert_v2ind_PoR(para_std_PoR,zero_index);

[PHIQ,muQ,delta1,delta0] = JSZ_transform(lamQ_ols_r,rinfQ_ols_r,sig_ols_r,mature,W,iota);
[B,A] = ABdiff(mature,delta1,PHIQ,delta0,muQ,sig_ols_r);
Yhat = repmat(A',T,1) + Y1*B';
rmse = sqrt(sum(sum((yields - Yhat).^2))/T/J)*1200;
disp(['RMSE = ', num2str(rmse),' pts']);

%% SECOND STAGE -- RESTRICTED unbiased
flag_mean = true;
[Phi_tilde, mu_tilde, V_tilde] = est_unb_var(Y1, 1, flag_mean, 5000, 1000, 50, 1, 50000);

mu_mu_ur = mu_tilde;
PHI_mu_ur = Phi_tilde;
Sigma_mu = chol(V_tilde)';
parameters = convert_ind2v_PoR(lamQ_ols_ur,rinfQ_ols_ur,sig_ols_ur, sig_lam0_ols_ur,sig_lam1_ols_ur,zero_index);

disp('optimization -- MU R');
parameters =  fminunc('MCS_JSZ_pc_rest',parameters ,options,zero_index, mature,W,mu_mu_ur,PHI_mu_ur,WA,WB,Sigma_mu,H1,H1v,H2,W2);

[lamQ_mu_r,rinfQ_mu_r,sig_mu_r, sig_lam0_mu_r,sig_lam1_mu_r] = convert_v2ind_PoR(parameters,zero_index);
[PHIQ_mu_r,muQ_mu_r] = JSZ_transform(lamQ_mu_r,rinfQ_mu_r,sig_mu_r,mature,W,iota);

mu_mu_r = sig_lam0_mu_r + muQ_mu_r;
PHI_mu_r = sig_lam1_mu_r + PHIQ_mu_r;

%-----------asymptotic standard errors for unbiased estimation-------------
[H1_mu, H1v_mu] = information_matrix(Y1(1:end-1,:),V_tilde,1);
std_P_para_mu = reshape(sqrt(diag(inv(H1_mu))),4,3)';
R_PoR = blkdiag(H1_mu,H1v_mu,H2);
GAMMA_PoR = numgrad('para_str2red_vec_PoR_pc',parameters,mature,W,W2,zero_index);
variance_PoR = inv(GAMMA_PoR'*R_PoR*GAMMA_PoR);
para_std_PoR = sqrt(diag(variance_PoR));
[lamQ_mu_r_std,rinfQ_mu_r_std,sig_mu_r_std, sig_lam0_mu_r_std,sig_lam1_mu_r_std] = convert_v2ind_PoR(para_std_PoR,zero_index);

[PHIQ,muQ,delta1,delta0] = JSZ_transform(lamQ_mu_r,rinfQ_mu_r,sig_mu_r,mature,W,iota);
[B,A] = ABdiff(mature,delta1,PHIQ,delta0,muQ,sig_mu_r);
Yhat = repmat(A',T,1) + Y1*B';
rmse = sqrt(sum(sum((yields - Yhat).^2))/T/J)*1200;
disp(['RMSE = ', num2str(rmse),' pts']);

%% SECOND STAGE - UNRESTRICTED UNBIASED

parameters = convert_ind2v(lamQ_ols_ur,rinfQ_ols_ur,Sigma_mu);

disp('optimization -- MU UR');
parameters =  fminunc('MCS_JSZ_pc',parameters ,options,mature,W,WA,WB,Sigma_mu,H1v,H2,W2);

[lamQ_mu_ur,rQ_mu_ur,sig_mu_ur] = convert_v2ind(parameters);
[PHIQ_mu_ur,muQ_mu_ur] = JSZ_transform(lamQ_mu_ur,rQ_mu_ur,sig_mu_ur,mature,W,iota);
mu_mu_ur = mu_tilde;
PHI_mu_ur = Phi_tilde;

%% calculate forward rates
% 5-to-10 year forward rates
[f_ols_ur, frn_ols_ur,ftp_ols_ur] = tp_fb(mu_ols_ur, PHI_ols_ur,yields,mature, W,lamQ_ols_ur, rinfQ_ols_ur, sig_ols_ur,5*12+1,10*12,iota);
[f_ols_r,frn_ols_r,ftp_ols_r] = tp_fb(mu_ols_r, PHI_ols_r, yields, mature, W,lamQ_ols_r, rinfQ_ols_r, sig_ols_r,5*12+1,10*12,iota);
[f_mu_r, frn_mu_r,ftp_mu_r] = tp_fb(mu_mu_r, PHI_mu_r, yields, mature, W, lamQ_mu_r, rinfQ_mu_r, sig_mu_r,5*12+1,10*12,iota);
[f_mu_ur, frn_mu_ur,ftp_mu_ur] = tp_fb(mu_mu_ur, PHI_mu_ur,yields,mature, W,lamQ_mu_ur, rQ_mu_ur, sig_mu_ur,5*12+1,10*12,iota);

%% Decomposition
% A = [f_ols_ur, frn_ols_ur,ftp_ols_ur, f_ols_r,frn_ols_r,ftp_ols_r, f_mu_r, frn_mu_r,ftp_mu_r, f_mu_ur, frn_mu_ur,ftp_mu_ur]*1200;

%% results
disp('*** OLS-UR***');
disp('Table 5: 1st-3rd Columns');
sig_lam0_ols_ur'*1200
sig_lam0_ols_ur_std'*1200
sig_lam1_ols_ur
sig_lam1_ols_ur_std
rinfQ_ols_ur*1200
rinfQ_ols_ur_std*1200
lamQ_ols_ur'
lamQ_ols_ur_std'
sig_ols_ur*1200
sig_ols_ur_std*1200

disp('*** OLS-R***');
disp('Table 5: 4th-6th Columns');
sig_lam0_ols_r'*1200
sig_lam0_ols_r_std'*1200
sig_lam1_ols_r
sig_lam1_ols_r_std
rinfQ_ols_r*1200
rinfQ_ols_r_std*1200
lamQ_ols_r'
lamQ_ols_r_std'
sig_ols_r*1200
sig_ols_r_std*1200

disp('*** BC-R***');
disp('Table 5: 7th-9th Columns');
sig_lam0_mu_r'*1200
sig_lam0_mu_r_std'*1200
sig_lam1_mu_r
sig_lam1_mu_r_std
rinfQ_mu_r*1200
rinfQ_mu_r_std*1200
lamQ_mu_r'
lamQ_mu_r_std'
sig_mu_r*1200
sig_mu_r_std*1200


disp('Table 6: Output');
% largest eigenvalues
fprintf('$max(eig(\\Phi))$  ');
fprintf(' & %6.4f ', [max(abs(eig(PHI_ols_ur))), max(abs(eig(PHI_ols_r))), max(abs(eig(PHI_mu_r)))]);
fprintf('\\\\ \n');

% half-life
max_irf = 12*20; % twenty years
irf_ols_ur = irf_var1(PHI_ols_ur, max_irf);
hl_ols_ur = find(irf_ols_ur>.5,1,'last');
irf_ols_r = irf_var1(PHI_ols_r, max_irf);
hl_ols_r = find(irf_ols_r>.5,1,'last');
irf_mu_r = irf_var1(PHI_mu_r, max_irf);
hl_mu_r = find(irf_mu_r>.5,1,'last');
fprintf('half-life        ');
fprintf(' &   %u ', [hl_ols_ur, hl_ols_r, hl_mu_r]);
fprintf('\\\\ \n');

% IRF at 5y
fprintf('IRF at 5y');
fprintf(' & %6.2f ', [irf_ols_ur(60), irf_ols_r(60), irf_mu_r(60)]);
fprintf('\\\\ \n');

% volatilities
fprintf('\\sigma(f_t^{61,120})');
fprintf(' & %4.3f ', 1200*[std(f_ols_ur), std(f_ols_r), std(f_mu_r)]);
fprintf('\\\\ \n');

fprintf('\\sigma(\\tilde{f}_t^{61,120}) ');
fprintf(' & %4.3f ', 1200*[std(frn_ols_ur), std(frn_ols_r), std(frn_mu_r)]);
fprintf('\\\\ \n');

fprintf('\\sigma(ftp_t^{61,120}) ');
fprintf(' & %4.3f ', 1200*[std(ftp_ols_ur), std(ftp_ols_r), std(ftp_mu_r)]);
fprintf('\\\\ \n');


%% plot

[year month day] = datevec(dates);
yearIndex = [1; year(2:length(year))~=year(1:length(year)-1)]==1;
yearIndex = yearIndex.*(mod(year,5)==0);

figure(1);
set(gcf, 'Units', 'inches')
set(gcf, 'Position', [1, 1, 6, 8])

min_yrange = -1;
max_yrange = ceil(1200*max([forw;frn_ols_ur;frn_ols_r;frn_mu_r;ftp_ols_ur;ftp_ols_r;ftp_mu_r]));

% recession data
startRec = datenum({'31-Jul-1990','31-Mar-2001','31-Dec-2007'});
endRec = datenum({'31-Mar-1991','30-Nov-2001','30-Jun-2009'});
nRec = 3;
y=[min_yrange max_yrange max_yrange min_yrange];
colorRec = [.7,.7,.7];

% plot risk-neutral rates
subplot('Position', [.08, .55, .9, .40]);
plot(dates, 1200*forw, 'k-', 'Linewidth',2);
hold on;
for i=1:nRec
    x=[startRec(i) startRec(i) endRec(i) endRec(i)];
    h_fill = fill(x,y,colorRec);
    set(h_fill,'EdgeColor',colorRec);
end
set(gca, 'Layer', 'top')
h1 = plot(dates, 1200*forw, 'k-', 'Linewidth',2);
h2 = plot(dates, 1200*frn_ols_ur, 'c-', 'Linewidth',2);
h3 = plot(dates, 1200*frn_ols_r, 'r:', 'Linewidth',2);
h4 = plot(dates, 1200*frn_mu_r, 'b--', 'Linewidth',2);
ylabel('Percent');
axis([min(dates) max(dates) min_yrange max_yrange]);
line([min(dates) max(dates)],[0 0],'Color','k');
title('Risk-neutral rates');
legend([h1, h2, h3, h4], {'forward rate', 'OLS-UR', 'OLS-R', 'BC-R'});
set(gca, 'XTick', dates(yearIndex==1));
set(gca, 'XTickLabel', year(yearIndex==1));

% plot term premium
subplot('Position', [.08, .05, .9, .40]);
hold on;
for i=1:nRec
    x=[startRec(i) startRec(i) endRec(i) endRec(i)];
    h_fill = fill(x,y,colorRec);
    set(h_fill,'EdgeColor',colorRec);
end
set(gca, 'Layer', 'top')
h1 = plot(dates, 1200*forw, 'k-', 'Linewidth',2);
h2 = plot(dates, 1200*ftp_ols_ur, 'c-', 'Linewidth',2);
h3 = plot(dates, 1200*ftp_ols_r, 'r:', 'Linewidth',2);
h4 = plot(dates, 1200*ftp_mu_r, 'b--', 'Linewidth',2);
ylabel('Percent');
axis([min(dates) max(dates) min_yrange max_yrange]);
line([min(dates) max(dates)],[0 0],'Color','k');
title('Forward term premia');
set(gca, 'XTick', dates(yearIndex==1));
set(gca, 'XTickLabel', year(yearIndex==1));
set(gcf, 'Color', 'w');

%% save decomposition (plot data) 
% [year,month] = datevec(dates);
% A = [f_ols_ur, frn_ols_ur,ftp_ols_ur, f_ols_r,frn_ols_r,ftp_ols_r, f_mu_r, frn_mu_r,ftp_mu_r, f_mu_ur, frn_mu_ur,ftp_mu_ur]*1200;
% csvwrite('S:\RiddellS\MBauer\BRW\brw_table5\figure2_plotdata.csv',A);
