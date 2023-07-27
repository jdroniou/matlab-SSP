function timestepping(zcase,ucase,bmm,dW,ncell,nvert,vertex,area,cg,cell_v,X,Mass,A,B_indices,I_indices,imesh,Ndt,dt,bdry_vert,mesh,h)
%global ucase zcase
tol=1e-8;
itermax=1000;
forcediff=1;
usol_idt=zeros(nvert,Ndt+1);
usol_idt(:,1)=X;
fid = fopen('results.txt','w');
% Relaxation parameter (1 for no relaxation)
relax = 1;
ITER=0;
num_updates=0;
Res=0;
zb=zeros(nvert,1);



for idt=1:Ndt;
    %str = sprintf('idt=%d / %d\n',idt,Ndt);
    % forkprint(fid,str);
   
    % Solution: non-linear iterations
    Xprev = X;

     %Source: Compute using (u_t)-laplace(zeta) and the exact solution
     b=assemble_source_stc(idt,cell_v,ncell,nvert,area,dW,Xprev,zcase);
     % Following is the source function for the deterministic SP
     % constructed from the exact solution.
%     b=assemble_source(idt*dt,cell_v,ncell,nvert,area,cg,zcase,ucase);
    b=b(I_indices);
    res_prev = 1e7;

    iter = 0;
    res = 1;

    rhs = Mass*Xprev + dt*b;
    %          ex_sol1=exact_solution(dt*idt,vertex)';
    %          myrhs=Mass*(exact_solution(dt*(idt-1),vertex)') + dt*b;
    %          check = norm(Mass(:,I_indices)*ex_sol1(I_indices)+forcediff*dt*A*zetau(ex_sol1)-myrhs, 2)
    %Dirichlet non-homogeneous BC
    zb(B_indices)=zetau((exact_solution(idt*dt,bdry_vert(B_indices,:),ucase))',zcase);%Zb(:,idt);%


    while (iter < itermax && res > tol)
        DeltaX=zeros(nvert,1);
        %%% Newton
        Nlin_V=ones(nvert,1);
        Nlin_V(I_indices)=dzeta(Xprev(I_indices),zcase);
        Nlin = spdiags(Nlin_V,0,nvert,nvert);
        Aglob = Mass +  forcediff*dt*A*Nlin;
        D1=Aglob(:,I_indices);
        D2=Aglob(:,B_indices);
        zu=Xprev;
        zu(I_indices)=zetau(zu(I_indices),zcase);
        nlsource = Mass * Xprev + forcediff*dt*A*zu;
        C2=zb(B_indices)-Xprev(B_indices);
        Rhs=rhs-nlsource-D2*C2;


        % bicgstab
        %             [L,U] = ilu(D1,struct('type','ilutp','droptol',1e-6));
        %             [deltaX,flag]=bicgstab(D1,Rhs,1e-6,20,L,U);
        %             if (flag ~= 0)
        %                 flag
        %                 error('bicgstab did not converge')
        %             end
        deltaX=D1\Rhs;

        DeltaX(I_indices)=deltaX;
        DeltaX(B_indices)=C2;

        X = Xprev+ relax*DeltaX;

        iter = iter+1;
        %residual by increments
        %res = norm(X-Xprev,Inf) / norm(Xprev,Inf);
        
        zetaX=X;
        zetaX(I_indices)=zetau(zetaX(I_indices),zcase);
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
        iter
        error('no convergence')
    end;
    usol = X;
    usol_idt(:,idt+1)=X;
    ITER=ITER+iter;
    Res=Res+abs(res);
    % ex_sol=exact_solution(idt*dt,vertex)';
    %Res = norm(ex_sol(B_indices)-X(B_indices),Inf) / norm(ex_sol(B_indices),Inf)

    % Write the solution and grid vtk files, to be plotted by "paraview"
    %         write_solution_vtk(usol,strcat('VTKout/solution',num2str(idt)),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);

    % Exact solution
%     ex_sol=exact_solution(dt*idt,vertex,ucase)';
    %         write_solution_vtk(ex_sol,strcat('VTKout/ex_sol',num2str(idt)),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
    %str = sprintf('Solution computed, iter=%d, res=%4.2e, max sol=%f, max ex_sol=%f\n', iter, res, max(usol(I_indices)), max(ex_sol(I_indices)));
    %forkprint(fid,str);
    %z_error=norm(zetau(usol(B_indices))-zetau(ex_sol(B_indices)),"inf")
end; % end time stepping
 ITER=ceil(ITER/Ndt);
 Res=Res/Ndt;
  num_updates=ceil(num_updates/Ndt);
  str = sprintf('Solution computed for bm=%d,num_time_steps=%d, Avg_iter=%d with Avg_relax=%d, Avg_res=%4.2e\n', bmm,Ndt, ITER,num_updates, Res);
  forkprint(fid,str);
%Saving solutions for each Brownian motion
save(strcat('solutions/BM',num2str(bmm),'msh',mesh(1:8),'tcuz',num2str(ucase),num2str(zcase)),'usol_idt','dt','Ndt','mesh','h');
end