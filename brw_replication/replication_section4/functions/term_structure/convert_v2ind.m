function [lamQ,rQ,sig_P] = convert_v2ind(para)
% [lamQ,rQ,sig_P] = convert_v2ind(para) converts a vector of parameters
% into individual parameters for JSZ's parameterization


para = para(:);
lamQ = para(1:3);
rQ = para(4);
sig_P = veclt(para(5:end))/100;