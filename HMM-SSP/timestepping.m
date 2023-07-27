function [usol,ex_sol]=timestepping(zc,uc,bmm,dW,aec,Area_CE,int_all,GBE_indices,BE_indices,all_cg,ncell,n,X,Mass,A,cell_e,Ndt,dt,mesh,h)
tol=1e-12;
itermax=1000;
forcediff=1;
usol_idt=zeros(n,Ndt+1);
usol_idt(:,1)=X;
fid = fopen('results.txt','w');
% Relaxation parameter (1 for no relaxation)
relax = 1;
ITER=0;
num_updates=0;
Res=0;
zb=zeros(n,1);



for idt=1:Ndt;
    % Solution: non-linear iterations
    Xprev = X;

    %Source (for deterministic SP): Compute using (u_t)-laplace(zeta) and the exact solution
    %b=assemble_source_dtr(idt*dt,Area_CE,all_cg,zc,cell_e,ncell,uc);

    %Source (for stochastic SP): Computing int(f(zeta(u))dW_t)
    b=assemble_source_stc(zc,idt,Area_CE,dW,Xprev);
    b=b(int_all);

    % str = sprintf('idt=%d / %d\n',idt,Ndt);
    % forkprint(fid,str);

    res_prev = 1e7;
    iter = 0;
    res = 1;

    rhs =Mass*Xprev+ dt*b;
    %         myrhs =Mass*(exact_solution(dt*(idt-1),all_cg)')+ dt*b;
    %         ex_sol12=exact_solution(dt*idt,all_cg)';
    %         check = norm(Mass*ex_sol12+forcediff*dt*A*zetau(ex_sol12,zc)-myrhs, 2)
    %Dirichlet non-homogeneous BC
    zb(GBE_indices)=zetau((exact_solution(idt*dt,aec(BE_indices,:),uc))',zc);

    DeltaX=zeros(n,1);
    while (iter < itermax && res > tol)
        %%% Newton
        Nlin_V=ones(n,1);
        Nlin_V(int_all)=dzeta(Xprev(int_all),zc);
        Nlin = spdiags(Nlin_V,0,n,n);
        Aglob = Mass +  forcediff*dt*A*Nlin;
        D1=Aglob(:,int_all);
        D2=Aglob(:,GBE_indices);
        zu=Xprev;
        zu(int_all)=zetau(zu(int_all),zc);
        nlsource = Mass * Xprev + forcediff*dt*A*zu;
        C2=zb(GBE_indices)-Xprev(GBE_indices);
        Rhs=rhs-nlsource-D2*C2;

        %bicgstab
%         [L,U] = ilu(D1,struct('type','ilutp','droptol',1e-6));
%         [deltaX,flag]=bicgstab(D1,Rhs,1e-10,20,L,U);
%         if (flag ~= 0)
%             flag
%             error('bicgstab did not converge')
%         end
        deltaX=D1\Rhs;

        DeltaX(int_all)=deltaX;
        DeltaX(GBE_indices)=C2;

        X = Xprev+ relax*DeltaX;

        iter = iter+1;
        % residual by increments
        %      res = norm(X-Xprev,Inf) / norm(Xprev,Inf);

        zetaX=X;
        zetaX(int_all)=zetau(zetaX(int_all),zc);
        res = norm(Mass*X+forcediff*dt*A*zetaX-rhs, 2);

        %% If the residual is too large, we don't update and we reduce relax
        if (res > 1.2*res_prev)
            iter = iter-1;
            relax = relax/5;
            num_updates=num_updates+1;
        else
            Xprev = X;

            res_prev = res;
            relax = min(1,relax*1.2);
        end;

    end; % end nonlinear iterations



    if (iter==itermax)
        res
        %iter
        error('no convergence')
    end;
    usol = X;
    usol_idt(:,idt+1)=X;
    ITER=ITER+iter;
    Res=Res+abs(res);

    % Write the solution and grid vtk files, to be plotted by "paraview"
    %         write_solution_vtk(usol,strcat('VTKout/solution',num2str(idt)),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);

    % Exact solution
    %         ex_sol=exact_solution(dt*idt,all_cg,uc)';
    %         write_solution_vtk(ex_sol,strcat('VTKout/ex_sol',num2str(idt)),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
    %str = sprintf('Solution computed, iter=%d, res=%4.2e, max sol=%f, max ex_sol=%f\n', iter, res, max(usol(int_all)), max(ex_sol(int_all)));
    %forkprint(fid,str);

end; % end time stepping

ITER=ceil(ITER/Ndt);
Res=Res/Ndt;
num_updates=ceil(num_updates/Ndt);
str = sprintf('Solution computed for bm=%d,num_time_steps=%d, Avg_iter=%d with Avg_relax=%d, Avg_res=%4.2e\n', bmm,Ndt, ITER,num_updates, Res);
forkprint(fid,str);
%Saving solutions for each Brownian motion
save(strcat('solutions/BM',num2str(bmm),'msh',mesh(1:8),'tcuz',num2str(uc),num2str(zc)),'usol_idt','dt','Ndt','mesh','h');
end