%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HMM code for u_t-Delta zeta(u)=f(zeta(u))*dW with Dirichlet BC
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
global t0;
global uc;
uc=1;
zcase=1;
forcediff=1;

t0=.1;


% Final time
T=1;

% Number of Brownian motions
nbm=2;

% Relaxation parameter (0 for no relaxation)
relax = 1;
% a paramenter to disturbute a portion of Mass from cells center to edges
r=0.5;
%%
% Sequence of meshes over which we want to run the scheme
% The meshes are available at https://github.com/jdroniou/HHO-Lapl-OM

%Hexagonal meshes
%meshes={'hexa1_01.mat';'hexa1_02.mat';'hexa1_03.mat'};%'hexa1_04.mat';'hexa1_05.mat'};

%Triangular meshes
meshes={'mesh1_01.mat';'mesh1_02.mat'};%'mesh1_03.mat';'mesh1_04.mat'};%'mesh1_05.mat'};
nbmeshes=size(meshes,1);

h=zeros(nbmeshes,1);
Ndt=zeros(nbmeshes,1);

% To see the results printed in file
fid = fopen('results.txt','w');

%% Start meshes
for imesh=1:nbmeshes
    % Load mesh
    loadmesh=strcat('load ../matlab_meshes/',meshes{imesh});
    str = sprintf('%s\n',loadmesh);
    forkprint(fid,str);
    eval(loadmesh);
    % Compute real centers of gravity and mesh diameter
    cg=gravity_centers(ncell,cell_v,vertex,area);
    h(imesh)=max(abs(diam))

    % Time steps

    Ndt(imesh)=ceil(T/h(imesh)^2)
    %     if (imesh>1)
    %       Ndt(imesh) = Ndt(imesh-1)*4;
    %     end;

    %% Distributing mass among edges
    area_edges=zeros(nedge,1);
    for i=1:ncell
        nbe=size(cell_e{i},2);
        rs=ones(nbe,1);

        area_edges(cell_e{i})=area_edges(cell_e{i})+((1-r)*area(i)).*rs/nbe;
    end;
    area_cells=r.*area;
    %% Initialise RHS and unknown
    ex_sol=zeros(ncell+nedge,1);
    ex_sol(1:ncell)=exact_solution(0,cg);
    for i=1:ncell
        nbe=size(cell_e{i},2);
        % midpoints of edges
        vertex_loc=vertex(cell_v{i},:);
        xs=zeros(nbe,2);
        xs([1:nbe],:)=(vertex_loc([1:nbe],:)+vertex_loc([2:nbe+1],:))/2;
        ex_sol(ncell+cell_e{i})=exact_solution(0,xs);
    end;

    X=ex_sol;%transform(ex_sol,ncell,'to_unknowns');

    %% ASSEMBLE MATRIX of Laplacian
    %
    [A,b]=assemble_diffusion_system(cell_v,cell_n,cell_e,ncell,nedge,vertex,area,center,cg);
    A=sparse(A);

    % Time steppings
    dt=T/Ndt(imesh);
    %% Loading brownian motion
    %%%%Constructing Brownian motion (Finacial toolbox is  required)
    W=zeros(Ndt(imesh),nbm);

    for bmm=1:nbm
        [Path,Time,dW]=simulate(bm(0,1),Ndt(imesh));
        W(:,bmm)=dW;
    end
    disp('Brownian motions loaded');
    %% Calculating BC
    fbe=zeros(nedge,1);
    xs=zeros(nedge,2);
    disp('Calculating BC and interior cells');

    % Dirichlet BC
    %  These BC are set up based on the fact that, on rows corresponding to boundary edges,
    %   A has 1 on the diagonal and 0 elsewhere. Since it's multiplied by forcediff * dt in the
    %   global system, so is the exact solution to define the BC

    for i=1:ncell
        I=find(cell_n{i}==0);
        if (size(I,2)>0)
            % midpoints of the boundary edges
            xs(cell_e{i}(I),:)=(vertex(cell_v{i}(I),:)+vertex(cell_v{i}(I+1),:))/2;
            fbe(cell_e{i}(I))=1;%% Fixing it in order to find the interior cells
        end;
    end
    % Calculating and saving RHS only for boundary edges along with their labels
    BE_indices=find(fbe);
    sizeBE=size(BE_indices,1);
    Rhs=zeros(sizeBE,1);
    for idt=1:Ndt(imesh);
        Rhs(:,idt)=area_edges(BE_indices).*exact_solution(idt*dt,xs(BE_indices,:))+forcediff*dt*zetau(exact_solution(idt*dt,xs(BE_indices,:)),zcase);
    end
    clear xs
    disp('Done')
    IE_indices=find(~fbe);
    z_norm=zeros(nbm,1);
    xi_sol=zeros(size(ncell+nedge,1),1);
    %% Timestepping for each Brownian motion using Parallel computing Toolbox
    % Computation is independend for each Brownian motion
    parfor bmm=1:nbm
        timestepping(bmm,zcase,Rhs,IE_indices,BE_indices,X,A,W(:,bmm),imesh,Ndt(imesh),ncell,nedge,dt,area_cells,area_edges,h,meshes{imesh});
    end
end; % end meshes

disp(['Elapsed time is ' num2str(round(toc/60,1)) ' minutes.'])
