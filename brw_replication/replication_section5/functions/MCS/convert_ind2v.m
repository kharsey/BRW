function [para] = convert_ind2v(lamQ,rinfQ,sig)
% [para] = convert_ind2v(lamQ,rQ,sig_P) converts JSZ's inidividual parameters
% into a vector.
%
% INPUT:    lamQ: diagonal mean reversion matrix
%           rQ: scalar, expectation of long term short rate under Q
%           sig_P: cholesky decomposition of covariance matrix
% Last updated 6/15/2010

para = [lamQ(:);rinfQ*10;ltvec(sig)*1000];