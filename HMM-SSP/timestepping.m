function timestepping(bmm,zc,Rhs,IE_indices,BE_indices,X,A,dW,imesh,Ndt,ncell,nedge,dt,area_cells,area_edges,h,mesh)
global uc;
n=ncell+nedge;
rhs=zeros(n,1);
zcase=zc;
fid = fopen('results.txt','a');
str = sprintf('Time step= %d,Brownian step= %d\n',Ndt,bmm);
forkprint(fid,str);

usol_idt=zeros(n,Ndt+1);
usol_idt(:,1)=X;

forcediff=1;
tol=1e-7;
itermax=500;

% Relaxation parameter (1 for no relaxation)
relax = 1;


%Counters
ITER=0;
num_updates=0;
max_iter=0;

for idt=1:Ndt;
  
    % Solution: non-linear iterations
    Xprev = X;
    %Source: Computing int(f(zeta(u))dW_t)
    b=assemble_source(zcase,idt,area_cells,area_edges,dW,Xprev);
    res_prev = 1e7;
    iter = 0;
    res = 1;
    rhs(ncell+BE_indices)=Rhs(:,idt);
    int_all=[1:ncell ncell+IE_indices']';
    rhs(int_all) = [area_cells;area_edges(IE_indices)].*Xprev(int_all)+ dt*b(int_all);

    rep=0;
while (iter < itermax && res > tol)
        %%% Newton
        Mass=spdiags([area_cells; area_edges],0,n,n);
        Nlin = spdiags(dzeta(Xprev,zcase),0,n,n);
        Aglob = Mass +  forcediff*dt*A*Nlin;
        nlsource = Mass * Xprev + forcediff*dt*A*zetau(Xprev,zcase);

        % bicgstab
        [L,U] = ilu(Aglob,struct('type','ilutp','droptol',1e-6));
        [deltaX,flag]=bicgstab(Aglob,rhs-nlsource,1e-6,20,L,U);
        if (flag ~= 0)
            flag
            error('bicgstab did not converge')
        end
        X = Xprev + relax*deltaX;

        iter = iter+1;
        % residual by increments
        %      res = norm(X-Xprev,Inf) / norm(Xprev,Inf);
        % residual of non-linear system
        %       if (m > 1)
        betau= zetau(X,zcase);
        res = norm( (Mass*X+forcediff*dt*A*betau)-rhs , 2);


        %% If the residual is too large, we don't update and we reduce relax
        if (res > 1.2*res_prev)
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
    max_iter=max(max_iter,iter);
    ITER=ITER+iter;

end; % end time stepping

num_updates=ceil(num_updates/Ndt);
ITER=ceil(ITER/Ndt);
str = sprintf('Solution computed for bm=%d,num_time_steps=%d, Avg_iter=%d,Max_iter=%d, with Avg_relax=%d, res=%4.2e\n', bmm,Ndt, ITER,max_iter,num_updates, res);
forkprint(fid,str);
%
save(strcat('solutions/BM',num2str(bmm),'msh',mesh(1:8),'tcuz',num2str(zc)),'usol_idt','dt','Ndt','mesh','h');

end
