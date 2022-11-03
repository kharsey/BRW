function [lamQ,rinfQ,sig, sig_lam0,sig_lam1] = convert_v2ind_PoR(para,zero_index)
% [lamQ,rQ,sig_P] = convert_v2ind(para) converts a vector of parameters
% into individual parameters for JSZ's parameterization

if nargin > 1
    para_ = nan(length(zero_index),1);
    para_(zero_index) = zeros(sum(zero_index),1);
    para_(isnan(para_)) = para;
else
    para_ = para;
end

[lamQ,rinfQ,sig] = convert_v2ind(para_(1:10));
sig_lam0 = para_(11:13)/1000;
sig_lam1 = reshape(para_(14:end),3,3);