function [KQ_1p,KQ_0p,rho_1p,rho_0p,penalty] = JSZ_transform(lamQ,rQ,sig_P,mature,W,iota)
%This function transforms the "fundamental" Q variables into the Q
%parameters under JSZ canonical rep with observable variables. The
%transformation follows Prop2 of JSZ(2009) on page of 8.
%Last update: 2/17/2011
                                  % for N = 3;

[Bx] = ABdiff(mature,iota,diag(lamQ));
Bx = W*Bx;
if rcond(Bx)<1e-15
    KQ_1p = [];
    KQ_0p = [];
    rho_1p= [];
    rho_0p=[];
    penalty = 1e50;
else
    penalty = 0;
    [DUMMY,Ax] = ABdiff(mature,iota, diag(lamQ),rQ,zeros(3,1), chol(Bx\sig_P*sig_P'/Bx')');   
    Ax = W*Ax;
    KQ_1p = Bx*(diag(lamQ)-eye(3))/Bx;
    KQ_0p = -KQ_1p*Ax;
    rho_1p = Bx'\iota;  
    rho_0p = rQ - Ax'*rho_1p;
    KQ_1p = KQ_1p + eye(3);
end
