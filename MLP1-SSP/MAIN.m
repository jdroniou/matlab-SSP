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
tic
clear;
format long;

forcediff=1;
% ucase is to select a function to compute IC and BC and zcase=\zeta(u)
ucase=1;
zcase=1;
% Number of Brownian motions
nbm=5;
% Final times
T=1;
%dt_initial=0.01;%0.001;

%% Select meshes
% Sequence of meshes over which we want to run the scheme
% These meshes are available at https://github.com/jdroniou/... TBC
meshes={'mesh1_01.mat';'mesh1_02.mat';'mesh1_03.mat'};%'mesh1_04.mat';'mesh1_05.mat';'mesh1_06.mat'};
nbmeshes=size(meshes,1);
h=zeros(nbmeshes,1);
Ndt=zeros(nbmeshes,1);
% To see the results printed in a file
fid = fopen('results.txt','w');
%% Start meshes
for imesh=1:nbmeshes

    % Load mesh
    loadmesh=strcat('load ../matlab_meshes/',meshes{imesh});
    str = sprintf('%s\n',loadmesh);
    forkprint(fid,str);
    eval(loadmesh);

    % Compute real centers of gravity, mesh diameter and area of dual mesh
    cg=gravity_centers(ncell,cell_v,vertex,area);
    h(imesh)=max(abs(diam));
    dualarea=compute_dualarea(area,ncell,nvert,cell_v);

    % Time steps
    Ndt(imesh)=ceil(T/h(imesh)^(2))
    %        if (imesh>1)
    %         Ndt(imesh) = Ndt(imesh-1)*4;
    %        end;
    %% Initialise RHS and unknown
    % Initial condition and exact solution
    ex_sol=exact_solution(0,vertex,ucase); % Exact solution at t0=0
    write_solution_vtk(ex_sol,strcat('VTKout/solution0'),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
    X=ex_sol;
    %% ASSEMBLE Diffusion MATRIX
    [A]=assemble_diffusion_system(cell_v,cell_n,ncell,nvert,vertex);
    % Time steppings
    dt=T/Ndt(imesh);
    %% Loading brownian motion
    %%%%Constructing Brownian motion (Finacial toolbox is  required)
    W=zeros(Ndt(imesh),nbm);
    for bmm=1:nbm
        [Path,Time,dW]=simulate(bm(0,1),Ndt(imesh));
        W(:,bmm)=dW;
    end
    disp('Brownian motions are loaded');

    %% Dirichlet BC
    fbc=zeros(size(vertex,1),1);
    bdry_vert=zeros(size(vertex,1),2);
    disp('Calculating solutions at the Boundary');
    for idt=1:Ndt(imesh)
        for i=1:ncell
            I=find(cell_n{i}==0);
            if (size(I,2)>0)
                bdry_vert_indices = [cell_v{i}(I) cell_v{i}(I+1)];
                bdry_vert(bdry_vert_indices,:) = vertex(bdry_vert_indices,:);
                fbc(bdry_vert_indices)=1;% Fixing it in order to find the interior cells
            end
        end
    end
    %Rhs
    disp('Done')
    % Calculating and saving RHS only for boundary cells along with their labels
    bdry_vert;
    IC_indices=find(~fbc);
    BC_indices=find(fbc);
    sizeBC_indices=size(BC_indices,1);
    Rhs=zeros(sizeBC_indices,1);
    for idt=1:Ndt(imesh)
        Rhs(:,idt) = dualarea(BC_indices).*(exact_solution(idt*dt,bdry_vert(BC_indices,:),ucase)) ...
            + forcediff*dt.*zetau((exact_solution(idt*dt,bdry_vert(BC_indices,:),ucase)),zcase);
    end
    %% Timestepping for each Brownian motion using Parallel computing Toolbox
    % Computation is independend for each Brownian motion
    z_norm=zeros(nbm,1);
    xi_sol=zeros(size(ncell+nedge,1),1);
    parfor bmm=1:nbm
        timestepping(bmm,vertex,ucase,zcase,Rhs,IC_indices,BC_indices,X,A,W(:,bmm),Ndt(imesh),cell_v,ncell,nvert,area,dualarea,dt,meshes{imesh},h);
    end
end% end meshes
disp(['Elapsed time is ' num2str(round(toc/60,1)) ' minutes for ' num2str(nbm) ' Brownian motions.' ])


   
