    function [X_sim] = genVAR(Phi, M, Y, p)
        % generate M data sets from VAR(p) model
        
        if nargin<4; p=1; end
        [T,k] = size(Y);
        
        Y_mean = mean(Y);
        Y_mean_rep = Y_mean' * ones(1,M);
        
        X_sim = zeros(k, T, M);
        
        if (p==1)
            % VAR(1)
            % 1. obtain residuals
            % use bootstrapped residuals
            resid = zeros(k, T-1);
            for t = 2:T
                resid(:,(t-1)) = Y(t,:)' - Y_mean' - (Phi * (Y(t-1,:) - Y_mean)' );
            end
            % 2. generate series
            % randomly select initial values from the data Y
            ind_start = randsample(T, M, true);
            X_sim(:,1,:) = Y( ind_start,:)';
            for t = 2:T
                u_sim = resid(:, randsample(T-1, M, true));
                X_sim(:,t,:) = Y_mean_rep + Phi * (squeeze(X_sim(:,t-1,:)) - Y_mean_rep) + u_sim;
            end
        elseif (p==2)
            % VAR(2)
            Phi1 = Phi(1:k,1:k);
            Phi2 = Phi(1:k, k+1:2*k);
            % randomly select initial values from the data Y
            ind = randsample(T-1, M, true);
            X_sim(:,1,:) = Y( ind,:)';
            X_sim(:,2,:) = Y( ind+1,:)';
            
            % use bootstrapped residuals
            % 1. obtain residuals
            resid = zeros(k, T-2);
            for t = 3:T
                resid(:,(t-2)) = Y(t,:)' - Y_mean' - (Phi1 * (Y(t-1,:) - Y_mean)' ) - (Phi2 * (Y(t-2,:) - Y_mean)' );
            end
            
            % 2. generate series
            for t = 3:T
                X_sim(:,t,:) = Y_mean_rep + Phi1 * (squeeze(X_sim(:,t-1,:)) - Y_mean_rep) + Phi2 * (squeeze(X_sim(:,t-2,:)) - Y_mean_rep) + resid(:, randsample(T-2, M, true));
            end
        else
            error('not implemented for p>2');
        end
    end