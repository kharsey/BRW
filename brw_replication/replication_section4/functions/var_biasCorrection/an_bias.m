function [bias] = an_bias(Phi, Sigma, T)
k = length(Phi);
% analytical bias approximation, based on Pope (1990)
Omega = Sigma * Sigma';
% long run variance Var(X_t)
Gamma0 = reshape(inv(eye(k^2) - kron(Phi, Phi))* Omega(:), k, k);
tmpsum = zeros(k,k);
lambda = eig(Phi');
for i=1:k
    tmpsum = tmpsum+lambda(i)*inv( eye(k) - lambda(i)*Phi' );
end
b = Omega * ( inv(eye(k)-Phi') + Phi'*inv(eye(k)-(Phi'*Phi')) + tmpsum ) * inv(Gamma0);

bias = -b/T;  % from Wang, Phillips, Yu
end