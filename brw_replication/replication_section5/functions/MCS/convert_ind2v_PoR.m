function [para] = convert_ind2v_PoR(lamQ,rinfQ,sig, sig_lam0,sig_lam1,zero_index)
% [para] = convert_ind2v(lamQ,rQ,sig_P) converts JSZ's inidividual parameters
% into a vector.
%
% INPUT:    lamQ: diagonal mean reversion matrix
%           rinfQ: scalar, expectation of long term short rate under Q
%           sig: cholesky decomposition of covariance matrix
%           sig_lam0,sig_lam1: prices of risk
%           zero_index: indexes indicating zero elements, which we need to
%                   take out in vector para
% Last updated 11/10/2011

para = [convert_ind2v(lamQ,rinfQ,sig);sig_lam0*1000;sig_lam1(:)];

if nargin > 5
    para(zero_index) = [];
end