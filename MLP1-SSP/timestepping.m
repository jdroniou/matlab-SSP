function timestepping(bmm,vertex,uc,zc,Rhs,IC_indices,BC_indices,X,A,dW,Ndt,cell_v,ncell,nvert,area,dualarea,dt,mesh,h)
ucase=uc;
zcase=zc;
fid = fopen('results.txt','a');
usol_idt=zeros(nvert,Ndt+1);
usol_idt(:,1)=X;
rhs=zeros(size(vertex,1),1);
forcediff=1;
% Tolerance
tol=1e-7;
itermax=100;
% Relaxation parameter (1 for no relaxation)
relax = 1;
%Counters
num_updates=0;
ITER=0;
for idt=1:Ndt
    %nonlinear iterations
    Xprev = X;
    %Source: Computing int(f(zeta(u))dW_t)
    b=assemble_source(zcase,idt,cell_v,ncell,nvert,area,dW,Xprev);
    res_prev = 1e7;
    iter = 0;
    res = 1;
    rhs(BC_indices)=Rhs(:,idt);
    rhs(IC_indices) = dualarea(IC_indices).*Xprev(IC_indices) + dt*b(IC_indices);

    while (iter < itermax && res > tol)
        %%% Newton
        Mass = spdiags(dualarea,0,nvert,nvert);%sparse(diag(dualarea));
        Nlin = spdiags(dzeta(Xprev,zcase),0,nvert,nvert);%sparse(diag(dzeta(Xprev,zcase)));
        Aglob = Mass +  forcediff*dt*A*Nlin;
        nlsource = Mass * Xprev + forcediff*dt*A*zetau(Xprev,zcase);

        % bicgstab
        [L,U] = ilu(Aglob,struct('type','ilutp','droptol',1e-6));
        [deltaX,flag]=bicgstab(Aglob,rhs-nlsource,1e-6,20,L,U);
        if (flag ~= 0)
            flag
            error('bicgstab did not converge')
        end
        %deltaX = Aglob\(rhs-nlsource);
        X = Xprev + relax*deltaX;


        iter = iter+1;
        % residual by increments
        %      res = norm(X-Xprev,Inf) / norm(Xprev,Inf);
        % residual of non-linear system
        res = norm((Mass*X+forcediff*dt*A*zetau(X,zcase))-rhs , 2);

        %% If the residual is too large, we don't update and we reduce relax
        if (res > 1.20*res_prev)

            iter = iter-1;
            relax = relax/5;
            num_updates=num_updates+1;
        else

            Xprev = X;
            res_prev = res;
            relax = min(1,relax*1.20);

        end;
     end; % end nonlinear iterations
    
    if (iter==itermax)
        res
        iter
        error('no convergence')
    end;
    usol_idt(:,idt+1)=X;
    ITER=ITER+iter;
end; % end time stepping
  ITER=ceil(ITER/Ndt);
  num_updates=ceil(num_updates/Ndt);
  str = sprintf('Solution computed for bm=%d,num_time_steps=%d, Avg_iter=%d with Avg_relax=%d, res=%4.2e\n', bmm,Ndt, ITER,num_updates, res);
  forkprint(fid,str);
 %Saving solutions for each Brownian motion 
 save(strcat('solutions/BM',num2str(bmm),'msh',mesh(1:8),'tcuz',num2str(ucase),num2str(zcase)),'usol_idt','dt','Ndt','mesh','h');
 %znorm=sum(Myznorm);
end