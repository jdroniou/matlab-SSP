%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HMM code for u_t-Delta u^m=0 with Dirichlet BC
%
%		Author: Jerome Droniou
%		Date: 04/01/20
%
% This code is provided as a companion of the article
%   "The gradient discretisation method for slow and fast diffusion porous media equations", J. Droniou and K. N. Le,
%   to appear in SIAM J. Numer. Anal. https://arxiv.org/abs/1905.01785
%
% Usage and modification of this code is permitted, but any scientific publication resulting from
% part of this code should cite the aforementioned article.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
format long;

% Parameters for nonlinearity and barenblatt solution
global t0;
global m;
global CB;
global uc;
uc=1;
zc=1;
forcediff=1;

t0=0.1;
%m=1;
CB=0.005;

%t0=.5;
%m=.3;
%CB=0.1;

% Final times
T=1;
% dt_initial=0.01;

% nonlinear iterations
tol=1e-7;
itermax=1000;
% Relaxation parameter (0 for no relaxation)
relax = 1;
r=0.5; % a paramenter to disturbute a portion of Mass from cells center to edges
%%
% Sequence of meshes over which we want to run the scheme
%%
% The meshes are available at https://github.com/jdroniou/HHO-Lapl-OM
meshes={'mesh1_1.mat';'mesh1_2.mat';'mesh1_3.mat'};%'mesh1_4.mat'};%'hexa1_3.mat';'hexa1_4.mat';'hexa1_5.mat'};

nbmeshes=size(meshes,1);
Lp1_error=zeros(nbmeshes,1);
L2zeta_error=zeros(nbmeshes,1);
h=zeros(nbmeshes,1);
Ndt=zeros(nbmeshes,1);

% To see the results printed in file
fid = fopen('results.txt','w');
str = sprintf('m=%f, t0=%f\n',m,t0);
forkprint(fid,str);
%%%fclose(fid);
% Ndt(1) = ceil(T/dt_initial);

for imesh=1:nbmeshes
    % Load mesh
    loadmesh=strcat('load ../../HHO-Lapl-OM-master/matlab_meshes/',meshes{imesh});
    str = sprintf('%s\n',loadmesh);
    forkprint(fid,str);
    eval(loadmesh);
    % Finding indices of boundary edges, interior edges
    % All edge centers
    aec=zeros(nedge,2);
    fbe=zeros(nedge,1);

    disp('Calculating BC and interior cells');
    for i=1:ncell
        % midpoints of all the boundary edges
        nbe=size(cell_e{i},2);
        aec(cell_e{i},:)=(vertex(cell_v{i}(1:nbe),:)+vertex(cell_v{i}(2:nbe+1),:))/2;
        I=find(cell_n{i}==0);
        if (size(I,2)>0)
            % Fixing it in order to find the boundary edges
            fbe(cell_e{i}(I))=1;
        end;
    end

    BE_indices=find(fbe);
    IE_indices=find(~fbe);
    int_all=[1:ncell ncell+IE_indices']';
    % Global boundary indices
    GBE_indices=ncell+BE_indices;
    n=ncell+nedge;

    % Compute real centers of gravity and mesh diameter
    cg=gravity_centers(ncell,cell_v,vertex,area);
    h(imesh)=max(abs(diam));
    % cell's gravity centers and all edge centers
    all_cg=[cg;aec];

    % Time steps
      Ndt(imesh)=ceil(T/h(imesh)^2);
%      if (imesh>1)
%          Ndt(imesh) = Ndt(imesh-1)*2;
%      end;

    str = sprintf('mesh= %s, h= %4.2e, time step= %4.2e \n',meshes{imesh},h(imesh),T/Ndt(imesh));
    forkprint(fid,str);
    %%Distributing mass among edges
    area_edges=zeros(nedge,1);
    for i=1:ncell
        nbe=size(cell_e{i},2);
        rs=ones(nbe,1);
        area_edges(cell_e{i})=area_edges(cell_e{i})+((1-r)*area(i)).*rs/nbe;
    end;
    area_cells=r.*area;
    Area_CE=[area_cells;area_edges];
    %% Initialise RHS and unknown
    ex_sol=zeros(n,1);
    ex_sol=exact_solution(0,all_cg,uc)';

    %     write_solution_vtk(ex_sol,strcat('VTKout/solution0'),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
    %     write_solution_vtk(ex_sol,'VTKout/ex_sol0',ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);

    X=ex_sol;%transform(ex_sol,ncell,'to_unknowns');

    X(GBE_indices)=zetau(ex_sol(GBE_indices),zc);
    %% ASSEMBLE MATRIX of Laplacian
    %
    [Diffmat,b]=assemble_diffusion_system(cell_v,cell_n,cell_e,ncell,nedge,vertex,area,center,cg);
    Mass = spdiags(Area_CE,0,n,n);
    Mass=Mass(int_all,:);
    A=Diffmat(int_all,:);
    % Time steppings
    dt=T/Ndt(imesh);
    ave_newton(imesh) = 0;
    zb=zeros(n,1);
    for idt=1:Ndt(imesh);

        %Source: Compute using (u_t)-laplace(zeta) and the exact solution
        b=assemble_source(idt*dt,Area_CE,all_cg,zc,cell_e,ncell,uc);
        b=b(int_all);

        % str = sprintf('idt=%d / %d\n',idt,Ndt(imesh));
        % forkprint(fid,str);
        % Solution: non-linear iterations
        Xprev = X;
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
                        [L,U] = ilu(D1,struct('type','ilutp','droptol',1e-6));
                        [deltaX,flag]=bicgstab(D1,Rhs,1e-6,20,L,U);
                        if (flag ~= 0)
                            flag
                            error('bicgstab did not converge')
                        end
%             deltaX=D1\Rhs;

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

            else
                Xprev = X;

                res_prev = res;
                relax = min(1,relax*1.2);
            end;

        end; % end nonlinear iterations


        ave_newton(imesh) = ave_newton(imesh) + iter;

        if (iter==itermax)
            res
            %iter
            error('no convergence')
        end;
        usol = X;


        % Write the solution and grid vtk files, to be plotted by "paraview"
        %         write_solution_vtk(usol,strcat('VTKout/solution',num2str(idt)),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);

        % Exact solution
        ex_sol=exact_solution(dt*idt,all_cg,uc)';
        %         write_solution_vtk(ex_sol,strcat('VTKout/ex_sol',num2str(idt)),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
        %str = sprintf('Solution computed, iter=%d, res=%4.2e, max sol=%f, max ex_sol=%f\n', iter, res, max(usol(int_all)), max(ex_sol(int_all)));
        %forkprint(fid,str);

    end; % end time stepping

    ave_newton(imesh) = ave_newton(imesh)/Ndt(imesh);

    % compute error
    [Lp1_error(imesh) L2zeta_error(imesh) H1zeta_error(imesh)] = compute_errors(usol,ex_sol,Area_CE,Diffmat,zc,int_all,GBE_indices);

    str = sprintf('Mesh %i. Errors: L^(2)=%4.2e, L^2 on zeta(u):%4.2e, H1 on zeta(u):%4.2e\n',imesh,Lp1_error(imesh),L2zeta_error(imesh),H1zeta_error(imesh),zc);
    forkprint(fid,str);

    % fraction of negative mass
    fraction_neg_mass(imesh) = abs( sum(area.*min(usol(1:ncell),0)) / sum(area.*abs(usol(1:ncell))) );
    str = sprintf('Mesh %i. Fraction negative mass %f\n',imesh,fraction_neg_mass(imesh));
    forkprint(fid,str);

    % nb of interior edges
    nedges_int(imesh) = size(IE_indices,1);


end; % end meshes


for imesh=1:nbmeshes-1 % convergence rate
    ocLp1(imesh)=log(Lp1_error(imesh)/Lp1_error(imesh+1))/log(h(imesh)/h(imesh+1));
    ocL2zeta(imesh)=log(L2zeta_error(imesh)/L2zeta_error(imesh+1))/log(h(imesh)/h(imesh+1));
    ocH1zeta(imesh)=log(H1zeta_error(imesh)/H1zeta_error(imesh+1))/log(h(imesh)/h(imesh+1));
end

str = sprintf('\nErrors in L^(2) and orders of convergence:\n');
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
fid = fopen('data_rates.dat','w');
fprintf(fid,'meshsize timestep L2error_zeta Lp1error H1error NedgesInt AveNewton FractionNegMass\n');
for i=1:nbmeshes
    fprintf(fid,'%f %f %f %f %f %d %f %f\n',h(i),T/Ndt(i),L2zeta_error(i),Lp1_error(i),H1zeta_error(i),nedges_int(i),ave_newton(i),fraction_neg_mass(i));
end;
fclose(fid);


