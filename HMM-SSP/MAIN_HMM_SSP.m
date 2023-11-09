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

% Parameters
tic
uc=1;
zc=1;
forcediff=1;

t0=0.1;

nbm=10;


% Final times
T=1;
% dt_initial=0.01;

r=0.5; % a paramenter to disturbute a portion of Mass from cells center to edges
%%
% Sequence of meshes over which we want to run the scheme
%%
% The meshes are available at https://github.com/jdroniou/HHO-Lapl-OM
%meshes={'hexa1_01.mat';'hexa1_02.mat';'hexa1_03.mat'};%'mesh1_4.mat'};%'hexa1_3.mat';'hexa1_4.mat';'hexa1_5.mat'};
meshes={'mesh1_01.mat';'mesh1_02.mat';'mesh1_03.mat';'mesh1_04.mat'};%'mesh1_05.mat';'mesh1_06.mat';'hexa1_5.mat'};

nbmeshes=size(meshes,1);
Lp1_error=zeros(nbmeshes,1);
L2zeta_error=zeros(nbmeshes,1);
h=zeros(nbmeshes,1);
Ndt=zeros(nbmeshes,1);

% To see the results printed in file
fid = fopen('results.txt','w');
str = sprintf('m=%f, t0=%f\n',2,t0);
forkprint(fid,str);
%%%fclose(fid);
% Ndt(1) = ceil(T/dt_initial);
%% Brownian motions for the finest mesh to be used for all other coarse meshes
%Reminder: we are using dt=h^2;
%Constructing Brownian motion (Finacial toolbox is  required)
fine_mesh=matfile(strcat('../matlab_meshes/',meshes{nbmeshes}));
h(nbmeshes)=max(abs(fine_mesh.diam));
Ndt(nbmeshes)=2*round(0.5*T/h(nbmeshes)^2);
dt=T/Ndt(nbmeshes);
WF=zeros(Ndt(nbmeshes),nbm);
     for bmm=1:nbm
         [Path,Time,dW]=simulate(bm(0,1),Ndt(nbmeshes));
         WF(:,bmm)=sqrt(dt)*dW;
     end
     save(strcat('BM/BM_mesh',num2str(nbmeshes)),'WF');
     clear WF;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for imesh=1:nbmeshes
    % Load mesh
    loadmesh=strcat('load ../matlab_meshes/',meshes{imesh});
    str = sprintf('%s\n',loadmesh);
    forkprint(fid,str);
    eval(loadmesh);
    % Finding indices of boundary edges, interior edges
    % All edge centers
    aec=zeros(nedge,2);
    fbe=zeros(nedge,1);

   
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
       Ndt(imesh)=2*round(0.5*T/h(imesh)^2);
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
        J=find(ismember(cell_e{i},BE_indices));
        I=size(J,2);
        if I>0
            rs(J)=0;
            nbe=nbe-I;
        end
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

    %% Loading brownian motion
    W=zeros(Ndt(imesh),nbm);
  load(strcat('BM/BM_mesh',num2str(nbmeshes)),'WF');
      if imesh<nbmeshes
        for i=1:Ndt(imesh)
            db=floor(Ndt(nbmeshes)/Ndt(imesh));
            W(i,:)=sum(WF([(db*(i-1)+1):(i*db)],:),1);

        end
        else
        W=WF;
     end
    clear WF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parfor bmm=1:nbm
        timestepping(zc,uc,bmm,W(:,bmm),aec,Area_CE,int_all,GBE_indices,BE_indices,all_cg,ncell,n,X,Mass,A,cell_e,Ndt(imesh),dt,meshes{imesh},h);
%     %Following error only works in the case of deterministic problem. Also, uncomment the correct source in the "timestepping". 
%     Sol=load(strcat('solutions/BM',num2str(bmm),'msh',meshes{imesh}(1:8),'tcuz',num2str(uc),num2str(zc)));;
%     Sol=Sol.usol_idt;
%     %Exact solution
%          ex_sol=exact_solution(dt*Ndt(imesh),all_cg,uc)';
%     usol=Sol(:,Ndt(imesh));
%     %compute error at the final times
%     [Lp1_error(imesh) L2zeta_error(imesh) H1zeta_error(imesh)] = compute_errors(usol,ex_sol,Area_CE,Diffmat,zc,int_all,GBE_indices);
%     str = sprintf('Mesh %i. Errors: L^(2)=%4.2e, L^2 on zeta(u):%4.2e, H1 on zeta(u):%4.2e\n',imesh,Lp1_error(imesh),L2zeta_error(imesh),H1zeta_error(imesh));
%     forkprint(fid,str);
    end
end; % end meshes
disp(['Elapsed time is ' num2str(toc) 'seconds' ])
