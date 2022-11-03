function [lamQ,rinfQ,sig] = convert_v2ind(para)
% [lamQ,rinfQ,sig] = convert_v2ind(para) converts a vector of parameters
% into individual parameters for JSZ's parameterization


para = para(:);
lamQ = para(1:3);
rinfQ = para(4)/10;
sig = veclt(para(5:end))/1000;