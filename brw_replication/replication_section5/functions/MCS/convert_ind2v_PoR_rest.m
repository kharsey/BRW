function [para] = convert_ind2v_PoR_rest(phiQ,rQ,sig, sig_lam0,sig_lam1,zero_index)
% [para] = convert_ind2v(lamQ,rQ,sig_P) converts JSZ's inidividual parameters
% into a vector.
%
% INPUT:    lamQ: diagonal mean reversion matrix
%           rQ: scalar, expectation of long term short rate under Q
%           sig_P: cholesky decomposition of covariance matrix
% Last updated 6/15/2010

%para = [convert_ind2v(phiQ,rQ,sig);sig_lam0(:)*1000;sig_lam1(:)];

% three elements zero - (3,1), (3,2), (2,2)
%para = [convert_ind2v(phiQ,rQ,sig);sig_lam0(:)*1000;sig_lam1(1:2,1);sig_lam1(1,2);sig_lam1(:,3)];

% only (2,2) element zero
%para = [convert_ind2v(phiQ,rQ,sig);sig_lam0(:)*1000;sig_lam1(:,1);sig_lam1(1,2);sig_lam1(3,2);sig_lam1(:,3)];

% (2,2) and (3,2)
%para = [convert_ind2v(phiQ,rQ,sig);sig_lam0(:)*1000;sig_lam1(:,1);sig_lam1(1,2);sig_lam1(:,3)];

% (3,2) and (2,3)
para = [convert_ind2v(phiQ,rQ,sig);sig_lam0(:)*1000;sig_lam1(:,1);sig_lam1(1:2,2);sig_lam1([1,3],3)];
