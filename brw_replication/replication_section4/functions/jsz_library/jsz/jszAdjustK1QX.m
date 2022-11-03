function [K1Q_X, isTypicalDiagonal] = jszAdjustK1QX(K1Q_X, eps0)
% function [K1Q_X, isTypicalDiagonal] = jszAdjustK1QX(K1Q_X, eps0);
%
% This function adjusts diagonal K1Q_X to give a non-diagonal but more
% computationally tractable K1Q_X.
%
%
% K1Q_X can fall into a few cases:
%   0. diagonal
%   1. not diagonal
%   2. zero eigenvalue
%   3. near repeated roots
% In cases 1-3, the diagonal closed form solver doesn't work, so compute differently.
%
% For case 2/3, we use the method in Joslin, Le, Singleton.  
%
% We order the diagonal of diagonal K1Q.
%
%
if nargin==1,
    eps0 = 1e-4;
end

diag_K1Q_X = diag(K1Q_X);
isDiagonal = all(all((K1Q_X==diag(diag_K1Q_X))));

if isDiagonal
    diag_K1Q_X = -sort(-diag_K1Q_X);
    K1Q_X = diag(diag_K1Q_X);
    
    hasNearUnitRoot = ~all(abs(diag_K1Q_X)>eps0); % Only applicable for diagonal
    hasNearRepeatedRoot = ~all(abs(diff(diag_K1Q_X))>eps0); % Only applicable for diagonal
    isTypicalDiagonal = isDiagonal & ~hasNearUnitRoot & ~hasNearRepeatedRoot;
else
    isTypicalDiagonal = false;
end


if isDiagonal && ~isTypicalDiagonal 
    inds = abs(diff(diag_K1Q_X))<eps0 | abs(diag_K1Q_X(1:end-1))<=eps0;
    if abs(diag_K1Q_X(end))<=eps0, 
        inds(end) = true;
    end
    K1Q_X(1:end-1,2:end) = K1Q_X(1:end-1,2:end) + diag(inds);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
