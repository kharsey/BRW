function [Phi_tilde, mu_tilde, V_tilde, Phi_sample] = est_unb_var(X, p, flag_mean, N, N_burn, B, check, B_check)
% est_unb_var - unbiased VAR estimation using stochastic approximation (SA)
% inputs:
%  X       REQUIRED. data matrix, T x k
%  p       order of the VAR, default VAR(1), p=2 possible, p>2 not implemented
%  flag_mean     flag whether mean- (TRUE) or median- (FALSE) unbiased estimation is desired
%           default TRUE (mean-unbiased)
%  N       number of iterations of the SA algorithm after burn-in (default 10,000)
%  N_burn  number of burn-in iterations (default 100)
%  B       number of bootstrap samples to calculate noisy measure of mean/median
%          of the OLS estimator (default 50)
%  check   flag whether closeness check is to be performed in the end (default TRUE)
%  B_check number of bootstrap samples for closeness check (default
%  100,000)
%
% outputs:
%  Phi_tilde  (companion form) estimated coefficient matrix
%  mu_tilde   estimated intercept
%  V_tilde    estimated variance-covariance matrix of residuals

% parameters
if nargin<8; B_check = 100000; end;
if nargin<7; check = true; end;
if nargin<6; B = 10; end;
if nargin<5; N_burn = 1000; end;
if nargin<4; N = 5000; end;
if nargin<3; flag_mean = true; end;
if nargin<2; p = 1; end;
gamma_i = .5;

if (flag_mean)
    fprintf('Mean-');
else
    fprintf('Median-');
end
fprintf('Unbiased VAR Estimation\n');
fprintf('N = %u, N_burn = %u, B = %u, B_check = %u, p = %u\n', [N,N_burn,B,B_check,p]);

[T,k] = size(X);

if (p==1)
    % first-order VAR
    
    % OLS
    [Phi_hat] = estVAR(X, 1, true, false);
    fprintf('largest absolute eigenvalue OLS:  %8.6f \n', max(abs(eig(Phi_hat))));
    
    % initialization for SA algorithm
    theta = zeros(k^2, N_burn+N);
    theta_hat = Phi_hat(:);
    theta(:,1) = theta_hat; % starting value
    
    % SA algorithm
    for j=1:N_burn+N-1
%        fprintf('****** iteration %u ******\n', j);
%        fprintf('eig '); fprintf('%8.6f ', eig(reshape(theta(:,j),k,k)));fprintf('\n');
        Phi_new = m_var(theta(:,j), B, 1, X, flag_mean);
%        fprintf('eig.new '); 
%        disp(eig(Phi_new));
        theta_new = Phi_new(:);
        d = theta_hat - theta_new;
        theta(:,j+1) = theta(:,j) + gamma_i*d;
        if ( (j>N_burn) && (mod(j-N_burn, 100) == 0) )
            % print some diagnostics
            fprintf('****** iteration %u ******\n', j-N_burn);
            theta_tilde = mean(theta(:,N_burn+1:j),2);
            Phi_tilde = reshape(theta_tilde,k,k);
            fprintf('largest absolute eigenvalue:  %8.6f \n', max(abs(eig(Phi_tilde))));
        end
    end
    
    theta_tilde = mean(theta(:,(N_burn+1):(N_burn+N)),2);
    
    if (check)
        % check whether mean/median of OLS is close to actual OLS estimates
        disp('... checking closeness of mean/median to actual estimates ...');
        [Phi_new, Phi_sample] = m_var(theta_tilde, B_check, 1, X, flag_mean, true);  % return sample
        theta_new = Phi_new(:);
        dist = sqrt(sum((theta_new - theta_hat).^2)/length(theta_new));
        fprintf('root mean square distance: %8.6f \n', dist);
    else
        Phi_sample = NaN;
    end
    
    % unbiased estimates
    Phi_tilde = reshape(theta_tilde,k,k);
    mu_tilde = (eye(k) - Phi_tilde) * mean(X)';
    
    % residuals and their variance-covariance matrix
    Xdem = X - ones(T,1)*mean(X);
    resid_tilde = Xdem(2:T,:)' - Phi_tilde * Xdem(1:(T-1),:)';
    V_tilde = resid_tilde * resid_tilde'/(T-1);
elseif (p==2)
    % second-order VAR
    
    % OLS
    [F_hat] = estVAR(X, 2, true, false);
    fprintf('largest absolute eigenvalue OLS:  %8.6f \n', max(abs(eig(F_hat))));
    
    % initialization for SA algorithm
    theta = zeros(2*k^2, N_burn+N);
    tmp = F_hat(1:k,:); theta_hat = tmp(:); 
    theta(:,1) = theta_hat; % starting value
    
    % SA algorithm
    for j=1:N_burn+N-1
        F_new = m_var(theta(:,j), B, 2, X, flag_mean);
        theta_new = F_new(:);
        d = theta_hat - theta_new;
        theta(:,j+1) = theta(:,j) + gamma_i*d;
        if ( (j>N_burn) && (mod(j-N_burn, 100) == 0) )
            % print some diagnostics
            fprintf('****** iteration %u ******\n', j-N_burn);
            theta_tilde = mean(theta(:,N_burn+1:j),2);
            F_tilde = [reshape(theta_tilde,k,2*k); [eye(k), zeros(k,k)]];
            fprintf('largest absolute eigenvalue:  %8.6f \n', max(abs(eig(F_tilde))));
        end
    end
    
    theta_tilde = mean(theta(:,(N_burn+1):(N_burn+N)),2);
    
    if (check)
        % check whether mean/median of OLS is close to actual OLS estimates
        disp('... checking closeness of mean/median to actual estimates ...');
        F_new = m_var(theta_tilde, B_check, 2, X, flag_mean);
        dist = sqrt(sum((F_new(:) - theta_hat).^2)/length(theta_new));
        fprintf('root mean square distance: %8.6f \n', dist);
    end
    
    % unbiased estimates
    Phi1_tilde = reshape(theta_tilde(1:k^2),k,k);
    Phi2_tilde = reshape(theta_tilde(k^2+1:2*k^2),k,k);
    F_tilde = [Phi1_tilde, Phi2_tilde; [eye(k), zeros(k,k)]];
    mu_tilde = (eye(k) - Phi1_tilde - Phi2_tilde) * mean(X)';
    Phi_tilde = F_tilde;
    
    % residuals and their variance-covariance matrix
    Xdem = X - ones(T,1)*mean(X);
    resid_tilde = Xdem(3:T,:)' - Phi1_tilde * Xdem(2:(T-1),:)' - Phi2_tilde * Xdem(3:T,:)';
    V_tilde = resid_tilde * resid_tilde'/(T-2);    
end

end
