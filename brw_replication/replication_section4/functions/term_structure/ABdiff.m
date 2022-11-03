function [B,A] = ABdiff(mature,delta1,rhoQ,delta0,cQ,sig) 
% ABdiff_Latent calculates A's and B's with difference equations
% USAGE: [A1,A2,B1,B2] = ABdiff_Latent(para_str, nN) 
% INPUTS:   para_str:   a structure of parameters
%           mature:     a structure containing mature.exact and
%                       mature.error
% OUTPUTS:  A1:         3*1 vector
%           A2:         Ne*1 vector
%           B1:         3*3 matrix
%           B2:         Ne*3 matrix
% last updated on April/5/10

nN = max(mature);

B_bar_temp = -delta1;
B = zeros(3,nN);
B(:,1) = -B_bar_temp;
if nargout <2
    for i = 2:nN
        B_bar_temp = rhoQ'* B_bar_temp - delta1;
        B(:,i) = -B_bar_temp/i;      
    end      
else
    A_bar_temp = -delta0;  %initial value for A_bar and B_bar
    A = zeros(nN,1);
    A(1)   = -A_bar_temp;  %A1 and B1
    for i = 2:nN
        A_bar_temp = A_bar_temp + B_bar_temp'*cQ + 1/2*(B_bar_temp'*sig*sig'*B_bar_temp)-delta0;
        B_bar_temp = rhoQ'* B_bar_temp - delta1;
        A(i)   = -A_bar_temp/i;
        B(:,i) = -B_bar_temp/i;      
    end 
    A = A(mature);
end
B = B';
B = B(mature,:);