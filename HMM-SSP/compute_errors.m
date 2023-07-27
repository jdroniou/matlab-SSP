% Compute the L2 error on the function u and the gradient of u^m
%
function [Lp1,L2zeta,H1beta]=compute_errors(Y,Yexact,weights,DiffMat,zc,int_all,GBE_indices)
Yexact(GBE_indices)=zetau(Yexact(GBE_indices),zc);
norm_p1_exact = sum(weights.*abs(Yexact).^(2))^(1/(2));
error_p1 = sum(weights.*abs(Yexact-Y).^(2))^(1/(2));
Lp1 = error_p1 / norm_p1_exact;

zetaY=Y;
zetaY(int_all) = zetau(Y(int_all),zc);
zetaYexact=Yexact;
zetaYexact(int_all) =  zetau(Yexact(int_all),zc);
norm_L2zeta_exact = sum(weights.*abs(zetau(Yexact,zc)).^(2))^(1/2);
error_L2zeta = sum(weights.*(abs(zetau(Yexact,zc))-abs(zetau(Y,zc))).^2 )^(1/2);
L2zeta = error_L2zeta / norm_L2zeta_exact;

% Energy

error_H1zeta = ((zetaY-zetaYexact)' * DiffMat * (zetaY-zetaYexact) )^(1/2);
norm_H1zeta_exact = ((zetaYexact)' * DiffMat * (zetaYexact) )^(1/2);
H1beta = error_H1zeta / norm_H1zeta_exact;

