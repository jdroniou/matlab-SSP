function f=fs(t,z,zcase,ucase)

% NOTE: there shouldn't be any "testcase" in the source term f1
% A proper way to compute f1 is to write, for f1=ut-(zeta'(u)*trace(Hessian of u)+zeta''(u)*norm(gradu)^2. The Hessian
% matrix of the exact solution is computed in exact_solution.
[u,ut,gradu,D2u]=exact_solution(t,z,ucase);
[Zu,gradZ,D2Z]=zetau(u,zcase);
f = ut-(trace(D2u)*gradZ+((norm(gradu))^2)*D2Z);
end