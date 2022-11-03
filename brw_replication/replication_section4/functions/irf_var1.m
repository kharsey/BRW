function irfvec = irf_var1(Phi, maxlag)
% calculate impulse response function for a VAR(1)
% for the first variable in response to shocks to the first variable

if nargin<2
    maxlag = 500;
end

if (length(Phi)>1)
    irfvec = zeros(maxlag, 1);
    Psi = Phi;
    irfvec(1) = Phi(1,1);
    for i=2:maxlag
        Psi = Phi * Psi;
        irfvec(i) = Psi(1,1);
    end
else
    error('scalar');
end
end