%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mass-lumped P1 code for u_t-Delta zeta(u)=f with Dirichlet BC
%
%    Author: Jerome Droniou
%    Date: 04/01/20
%
% This code is provided as a companion of the article
%   "The gradient discretisation method for slow and fast diffusion porous media equations", J. Droniou and K. N. Le,
%   to appear in SIAM J. Numer. Anal. https://arxiv.org/abs/1905.01785
%
%
% Usage and modification of this code is permitted, but any scientific publication resulting from
% part of this code should cite the aforementioned article.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

format long;
tic
% Parameters for nonlinearity and barenblatt solution
global t0;
global CB;
%global ucase zcase;
ucase=1;
zcase=1;

forcediff=1;

t0=0.1;
CB=0.005;
%t0=.5;
%CB=0.1;

% Final times
T=1;
dt_initial=0.1;%0.001;

% nonlinear iterations
tol=1e-8;
itermax=1000;
% Relaxation parameter (1 for no relaxation)
relax = 1;

%%
% Sequence of meshes over which we want to run the scheme
%%
% The meshes are available at https://github.com/jdroniou/HHO-Lapl-OM
meshes={'mesh1_1.mat';'mesh1_2.mat';'mesh1_3.mat'};%'mesh1_4.mat'};%'mesh1_5.mat'};%'mesh1_6.mat'};%'mesh1_7.mat'};

nbmeshes=size(meshes,1);
Lp1_error=zeros(nbmeshes,1);
L2zeta_error=zeros(nbmeshes,1);
h=zeros(nbmeshes,1);
Ndt=zeros(nbmeshes,1);

% To see the results printed in file
fid = fopen('results.txt','w');
str = sprintf('t0=%f\n',t0);%sprintf('m=%f, t0=%f\n',m,t0);
forkprint(fid,str);
%%%fclose(fid);
Ndt(1) = ceil(T/dt_initial);

for imesh=1:nbmeshes
    % Load mesh
    loadmesh=strcat('load ../../HHO-Lapl-OM-master/matlab_meshes/',meshes{imesh});
    str = sprintf('%s\n',loadmesh);
    forkprint(fid,str);
    eval(loadmesh);
    % Compute real centers of gravity, mesh diameter and area of dual mesh
    cg=gravity_centers(ncell,cell_v,vertex,area);
    h(imesh)=max(abs(diam));%
   
    %% Finding boundary vertices
    fbc=zeros(size(vertex,1),1);
    bdry_vert=zeros(size(vertex,1),2);
    for i=1:ncell
        I=find(cell_n{i}==0);
        if (size(I,2)>0)
            bdry_vert_indices = [cell_v{i}(I) cell_v{i}(I+1)];
            bdry_vert(bdry_vert_indices,:)  = vertex(bdry_vert_indices,:);
            fbc(bdry_vert_indices)=1;
        end
    end
    I_indices=find(~fbc);
    B_indices=find(fbc);
    nbvert=size(B_indices,1);

    dualarea=compute_dualarea(area,ncell,nvert,cell_v,B_indices);
    %Finding the associated cells with the boundary vertices
%     cellwithBI=zeros(nbvert,5);
%     bi=1;
%     while bi<=nbvert
%         Bi=B_indices(bi);
%         cellwithBI(bi,1)=Bi;
%         p=2;
%         for i=1:ncell
%             for j=1:3
%                 if cell_v{i}(j)==Bi
%                     cellwithBI(bi,p)=i;
%                     p=p+1;
%                 end
%             end
%         end
%         bi=bi+1;
%     end
% 
% cellwithBI


    % Distrubting area of boundary vertices among interior vertices
    % equally and puting Mass 0 at boundary vertices
%     dualarea=zeros(nvert,1);
%     dualarea(I_indices)=Dualarea(I_indices)+sum(Dualarea(B_indices))/sum(Dualarea(I_indices))*Dualarea(I_indices);
    


    % Time steps
    %Ndt(imesh)=ceil(T/h(imesh)^2);
    if (imesh>1)
        Ndt(imesh) = Ndt(imesh-1)*2;
    end;

    str = sprintf('mesh= %s, h= %4.2e, time step= %4.2e \n',meshes{imesh},h(imesh),T/Ndt(imesh));
    forkprint(fid,str);

    %% Initialise RHS and unknown
    % Initial condition and exact solution
    ex_sol=exact_solution(0,vertex,ucase)'; % Exact solution at t0=0

%     write_solution_vtk(ex_sol,strcat('VTKout/solution0'),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
%     write_solution_vtk(ex_sol,'VTKout/ex_sol0',ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
    % Unknows at the boundary are "zeta(u)" while "u" at the interior points
    X=ex_sol;
    X(B_indices)=zetau(ex_sol(B_indices),zcase);

    %% ASSEMBLE MATRIX of Laplacian
    [Diffmat,b]=assemble_diffusion_system(cell_v,ncell,nvert,vertex);
    Mass = spdiags(dualarea,0,nvert,nvert);
    A=Diffmat(I_indices,:);
    Mass=Mass(I_indices,:);


    % Time steppings
    dt=T/Ndt(imesh);
    ave_newton(imesh) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ITER=0;
    num_updates=0;
    Res=0;
    zb=zeros(nvert,1);
    DeltaX=zeros(nvert,1);
    for idt=1:Ndt(imesh);
        %str = sprintf('idt=%d / %d\n',idt,Ndt(imesh));
        % forkprint(fid,str);
        %Source: Compute using (u_t)-laplace(zeta) and the exact solution
        b=assemble_source(idt*dt,cell_v,ncell,nvert,area,cg,zcase,ucase);
        b=b(I_indices);
        %b=0;
        % Solution: non-linear iterations
        Xprev = X;
        res_prev = 1e7;

        iter = 0;
        res = 1;

        rhs = Mass*Xprev + dt*b;
%          ex_sol1=exact_solution(dt*idt,vertex)';
%          myrhs=Mass*(exact_solution(dt*(idt-1),vertex)') + dt*b;
%          check = norm(Mass(:,I_indices)*ex_sol1(I_indices)+forcediff*dt*A*zetau(ex_sol1)-myrhs, 2)
        %Dirichlet non-homogeneous BC
        zb(B_indices)=zetau((exact_solution(idt*dt,bdry_vert(B_indices,:),ucase))',zcase);

        
        while (iter < itermax && res > tol)
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
        ave_newton(imesh) = ave_newton(imesh) + iter;

        if (iter==itermax)
            res
            iter
            error('no convergence')
        end;
        usol = X;
        ITER=ITER+iter;
        Res=Res+abs(res);
       % ex_sol=exact_solution(idt*dt,vertex)';
        %Res = norm(ex_sol(B_indices)-X(B_indices),Inf) / norm(ex_sol(B_indices),Inf)

        % Write the solution and grid vtk files, to be plotted by "paraview"
%         write_solution_vtk(usol,strcat('VTKout/solution',num2str(idt)),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);

        % Exact solution
        ex_sol=exact_solution(dt*idt,vertex,ucase)';
%         write_solution_vtk(ex_sol,strcat('VTKout/ex_sol',num2str(idt)),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
        str = sprintf('Solution computed, iter=%d, res=%4.2e, max sol=%f, max ex_sol=%f\n', iter, res, max(usol(I_indices)), max(ex_sol(I_indices)));
         forkprint(fid,str);
        %z_error=norm(zetau(usol(B_indices))-zetau(ex_sol(B_indices)),"inf")
    end; % end time stepping

    ITER=ceil(ITER/Ndt(imesh));
    Res=(Res/Ndt(imesh));
    num_updates=ceil(num_updates/Ndt(imesh));
    str = sprintf('Mesh %i, num_time_steps=%d, Avg_iter=%d with Avg_relax=%d, Avg_residue=%4.2e\n',imesh,Ndt(imesh),ITER,num_updates,Res);
    forkprint(fid,str);
    clear ITER Res num_updates

    ave_newton(imesh) = ave_newton(imesh)/Ndt(imesh);
   
    % compute error
    [Lp1_error(imesh) L2zeta_error(imesh) H1zeta_error(imesh)] = compute_errors(usol,ex_sol,dualarea,Diffmat,I_indices,zcase,ucase);

    str = sprintf('Mesh %i. Errors: L^(2)=%4.2e, L^2 on zeta(u):%4.2e, H1 on zeta(u):%4.2e\n',imesh,Lp1_error(imesh),L2zeta_error(imesh),H1zeta_error(imesh));
    forkprint(fid,str);


    % fraction of negative mass
    fraction_neg_mass(imesh) = abs( sum(dualarea.*min(usol,0)) / sum(dualarea.*abs(usol)) );
    str = sprintf('Mesh %i. Fraction negative mass %f\n',imesh,fraction_neg_mass(imesh));
    forkprint(fid,str);

    % nb of interior vertices
    nvert_int(imesh) = nvert;
    for i=1:ncell
        for j=1:size(cell_e{i},2)
        	  if (cell_n{i}(j)==0)
                  nvert_int(imesh) = nvert_int(imesh)-1;
              end
        end
    end
    Time(imesh)=toc;
end; % end meshes
disp(['Elapsed time is ' num2str(toc) 'seconds' ])
% convergence rate
for imesh=1:nbmeshes-1
    ocLp1(imesh)=log(Lp1_error(imesh)/Lp1_error(imesh+1))/log(h(imesh)/h(imesh+1));
    ocL2zeta(imesh)=log(L2zeta_error(imesh)/L2zeta_error(imesh+1))/log(h(imesh)/h(imesh+1));
    ocH1zeta(imesh)=log(H1zeta_error(imesh)/H1zeta_error(imesh+1))/log(h(imesh)/h(imesh+1));
end

str = sprintf('Errors in L^(2) and orders of convergence:\n');
forkprint(fid,str);
for imesh=1:nbmeshes
    if (imesh==1)
        str = sprintf('\t%4.2e\n',Lp1_error(imesh));
        forkprint(fid,str);
    else
        str = sprintf('\t%4.2e \t %4.2e\n',Lp1_error(imesh),ocLp1(imesh-1));
        forkprint(fid,str);
    end
end

str = sprintf('\nErrors in L^2 on zeta(u) and orders of convergence:\n');
forkprint(fid,str);
for imesh=1:nbmeshes
    if (imesh==1)
        str = sprintf('\t%4.2e\n',L2zeta_error(imesh));
        forkprint(fid,str);
    else
        str = sprintf('\t%4.2e \t %4.2e\n',L2zeta_error(imesh),ocL2zeta(imesh-1));
        forkprint(fid,str);
    end
end

str = sprintf('\nErrors in H1 on zeta(u) and orders of convergence:\n');
forkprint(fid,str);
for imesh=1:nbmeshes
    if (imesh==1)
        str = sprintf('\t%4.2e\n',H1zeta_error(imesh));
        forkprint(fid,str);
    else
        str = sprintf('\t%4.2e \t %4.2e\n',H1zeta_error(imesh),ocH1zeta(imesh-1));
        forkprint(fid,str);
    end
end

fclose(fid);

% Write data file
fid = fopen('data_rates_SP.dat','w');
fprintf(fid,'meshsize timestep L2error_zeta Lmp1error H1error NvertInt Time\n');
for i=1:nbmeshes
    fprintf(fid,'%f %f %f %f %f %d %f\n',h(i),T/Ndt(i),L2zeta_error(i),Lp1_error(i),H1zeta_error(i),nvert_int(i),Time(i));
end;
fclose(fid);
