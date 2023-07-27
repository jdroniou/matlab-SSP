% Compute the L2 error on the function zeta(u) and the gradient of zeta(u)
%
function [Lp1,L2zeta,H1zeta]=compute_errors(Y,Yexact,weights,DiffMat,I_indices,zcase,ucase)


% norm_p1_exact = sum(weights.*abs(Yexact).^(2))^(1/(2));
% error_p1 = sum(weights.*abs(Yexact-Y).^(2))^(1/(2));
% Lp1 = error_p1 / norm_p1_exact;

norm_p1_exact = sum(weights(I_indices).*abs(Yexact(I_indices)).^(2))^(1/(2));
error_p1 = sum(weights(I_indices).*abs(Yexact(I_indices)-Y(I_indices)).^(2))^(1/(2));
Lp1 = error_p1 / norm_p1_exact;

zetaY=Y;
zetaY(I_indices) = zetau(Y(I_indices),zcase);

zetaYexact =  zetau(Yexact,zcase);

norm_L2zeta_exact = sum(weights.*abs(zetaYexact).^(2))^(1/2);
error_L2zeta = sum(weights.*(abs(zetaYexact)-abs(zetaY)).^2 )^(1/2);
L2zeta = error_L2zeta / norm_L2zeta_exact;

% Energy

error_H1zeta = ((zetaY-zetaYexact)' * DiffMat * (zetaY-zetaYexact) )^(1/2);
norm_H1zeta_exact = ((zetaYexact)' * DiffMat * (zetaYexact) )^(1/2);
H1zeta = error_H1zeta / norm_H1zeta_exact;

