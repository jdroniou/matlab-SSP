%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mass-lumped P1 code for u_t-Delta zeta(u)=f(zeta(u))*dW with Dirichlet BC
%
%    Author: Muhammad Awais Khan
%    Date: 14/06/2023
% % %
% This code is provided as a companion of the article
%   "Numerical analysis of the stochastic Stefan problem",
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


nbm=1;
t0=0.1;
CB=0.005;
%t0=.5;
%CB=0.1;

% Final times
T=1;
%dt_initial=0.1;%0.001;

% nonlinear iterations



%%
% Sequence of meshes over which we want to run the scheme
%%
% The meshes are available at https://github.com/jdroniou/HHO-Lapl-OM
meshes={'mesh1_01.mat';'mesh1_02.mat';'mesh1_03.mat'};%'mesh1_04.mat';'mesh1_05.mat'};%'mesh1_06.mat'};

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
%Ndt(1) = ceil(T/dt_initial);

for imesh=1:nbmeshes
    % Load mesh
    loadmesh=strcat('load ../matlab_meshes/',meshes{imesh});
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
    % Time steps
    Ndt(imesh)=ceil(T/h(imesh)^2);
    %     if (imesh>1)
    %         Ndt(imesh) = Ndt(imesh-1)*2;
    %     end;

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

    %     Zb=zeros(nbvert,Ndt(imesh));
    %     for idt=1:Ndt(imesh)
    %          %dualarea(BC_indices).*(exact_solution(idt*dt,bdry_vert(BC_indices,:),ucase)) ...
    %           %  + forcediff*dt.*zetau((exact_solution(idt*dt,bdry_vert(BC_indices,:),ucase)),zcase);
    %         Zb(:,idt)=zetau((exact_solution(idt*dt,bdry_vert(B_indices,:),ucase))',zcase);
    %     end


    %% Loading brownian motion
    %%%%Constructing Brownian motion (Finacial toolbox is  required)
    W=zeros(Ndt(imesh),nbm);
    for bmm=1:nbm
        [Path,Time,dW]=simulate(bm(0,1),Ndt(imesh));
        W(:,bmm)=dW;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parfor bmm=1:nbm
        timestepping(zcase,ucase,bmm,W(:,bmm),ncell,nvert,vertex,area,cg,cell_v,X,Mass,A,B_indices,I_indices,imesh,Ndt(imesh),dt,bdry_vert,meshes{imesh},h);
    % Following error only works in the case of deterministic problem. Also, uncomment the correct source in the "timestepping". 
%     Sol=load(strcat('solutions/BM',num2str(bmm),'msh',meshes{imesh}(1:8),'tcuz',num2str(ucase),num2str(zcase)));
%     Sol=Sol.usol_idt;
%     % Exact solution
%     ex_sol=exact_solution(dt*Ndt(imesh),vertex,ucase)';
%     usol=Sol(:,Ndt(imesh));
% 
%     % compute error at the final times
%     [Lp1_error(imesh) L2zeta_error(imesh) H1zeta_error(imesh)] = compute_errors(usol,ex_sol,dualarea,Diffmat,I_indices,zcase,ucase);
%     str = sprintf('Mesh %i. Errors: L^(2)=%4.2e, L^2 on zeta(u):%4.2e, H1 on zeta(u):%4.2e\n',imesh,Lp1_error(imesh),L2zeta_error(imesh),H1zeta_error(imesh));
%     forkprint(fid,str);
    end




end; % end meshes
disp(['Elapsed time is ' num2str(toc) 'seconds' ])




