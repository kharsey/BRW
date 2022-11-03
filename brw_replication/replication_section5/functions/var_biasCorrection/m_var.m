function [Phi_new, Phi_sample] = m_var(theta, M, p, Y, flag_mean, flag_sample)
% find mean/median of OLS when DGP is VAR(p)
% works for VAR(1) and VAR(2) only
% inputs:
%  theta - vector of parameters determining Phi_1 (and Phi_2)
%  M - number of Monte Carlo replications
%  p - order of the VAR
%  Y - data, T x k matrix
%  flag_mean - flag whether mean (TRUE) or median (FALSE) is to be returned

[T, k] = size(Y);

if nargin<6; flag_sample=false; end;

if flag_sample
    Phi_sample = zeros(k,k,M);
else
    Phi_sample = NaN;
end;

if (p==1)
    % get coefficient matrix
    Phi_tilde =  reshape(theta, k, k);
    
    % simulate M datasets
    X_sim = genVAR(Phi_tilde, M, Y, p);
    
    % estimate VAR(1) on each series
    theta_new_i = zeros(M, k^2);
    for m=1:M
        Phi_new = estVAR(X_sim(:,:,m)', 1, true, false);
        if (flag_sample)
            Phi_sample(:,:,m) = Phi_new;
        end;
        theta_new_i(m,:) = Phi_new(:);
    end
    
    if (flag_mean)
        % find mean of OLS estimates
        Phi_new = reshape( mean(theta_new_i), k, k );
    else
        % find median of OLS estimates
        Phi_new = reshape( median(theta_new_i), k, k );
    end
    
elseif (p==2)
    % get coefficient matrix
    Phi_tilde =  reshape(theta, k, 2*k);
    
    % simulate M datasets
    X_sim = genVAR(Phi_tilde, M, Y, p);
    
    % estimate VAR(2) on each series
    theta_new_i = zeros(M, 2*k^2);
    for m=1:M
        F_new = estVAR(X_sim(:,:,m)', 2, true, false);
        tmp = F_new(1:k,:);
        theta_new_i(m,:) = tmp(:);
    end
    
    if (flag_mean)
        % find mean of OLS estimates
        Phi_new = reshape( mean(theta_new_i), k, 2*k );
    else
        % find median of OLS estimates
        Phi_new = reshape( median(theta_new_i), k, 2*k );
    end
else
    error('not yet implemented for order>2');
end

%        fprintf('Phi.new ');fprintf('%8.6f ', Phi_new); fprintf('\n');
%        fprintf('eig.new ');fprintf('%8.6f ', eig(Phi_new)); fprintf('\n');

end
