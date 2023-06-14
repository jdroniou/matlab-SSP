%In the following code, we will be calculating error in relavant norms for
% the reconstructed functions of \zeta(u), \nabla\zeta(u) and \Xi(u)

clear
format long;
tic
ucase=1;
zcase=1;
T=1;
nbm=5;

%% Select meshes
meshes={'mesh1_01.mat';'mesh1_02.mat';'mesh1_03.mat'};%'mesh1_04.mat';'mesh1_05.mat';'mesh1_06.mat'};
nbmeshes=size(meshes,1);
% F stands for the finer mesh
F=load(strcat('../matlab_meshes/',meshes{nbmeshes}));
%To compute diffusion matrix for norms, we do not regard boundary conditions! 
F_A=assemble_diffusion_system_forNorms(F.cell_v,F.ncell,F.nvert,F.vertex);
F_h=max(abs(F.diam));
F_Ndt=ceil(T/F_h^2);
F_dt=T/F_Ndt;
F_dualarea=compute_dualarea(F.area,F.ncell,F.nvert,F.cell_v);
%
Fst_znorm=zeros(nbm,1);
GFst_znorm=zeros(nbm,1);
Fst_Xi=zeros(nbm,1);
Fznorm=zeros(F_Ndt+1,1);
Est_eznorm=zeros(nbmeshes-1,1);
GEst_eznorm=zeros(nbmeshes-1,1);
Est_eXi=zeros(nbmeshes-1,1);
L2zeta=zeros(nbmeshes,1);
H1zeta=zeros(nbmeshes,1);
L1Xi=zeros(nbmeshes,1);
C_dof=zeros(nbmeshes-1,1);
C_h=zeros(nbmeshes-1,1);

fid = fopen('results.txt','w');

for imesh=1:nbmeshes-1
    
    str = sprintf('Loading %s and %s\n',meshes{imesh},meshes{nbmeshes});
    forkprint(fid,str);

    % C is for coarse mesh
    C=load(strcat('../matlab_meshes/',meshes{imesh}));
    C_A=assemble_diffusion_system_forNorms(C.cell_v,C.ncell,C.nvert,C.vertex);
    C_h(imesh)=max(abs(C.diam));
    C_Ndt=ceil(T/C_h(imesh)^2);
    C_dt=T/C_Ndt;
    C_dof(imesh)=C.nvert;
   
    Gst_eznorm=zeros(nbm,1);
    st_eznorm=zeros(nbm,1);
    st_eXi=zeros(nbm,1);
    Cst_Xi=zeros(nbm,1);
    ICTS_T=zeros(nbm,1);
    err_znorm=zeros(F_Ndt+1,1);
    Gerr_znorm=zeros(F_Ndt+1,1);
    Gerr_znorm2=zeros(F_Ndt+1,1);
%Computing for each Brownian motion       
    for bmm=1:nbm
        str = sprintf('mesh_0%d, b=%d / %d\n',C_Ndt,bmm,nbm);
        forkprint(fid,str);
        %Loading solutions of the coarsive mesh and the solutions on the
        %finest possible mesh as a reference solution
        CS=load(strcat('solutions/BM',num2str(bmm),'msh',meshes{imesh}(1:8),'tcuz',num2str(ucase),num2str(zcase)));
        FS=load(strcat('solutions/BM',num2str(bmm),'msh',meshes{nbmeshes}(1:8),'tcuz',num2str(ucase),num2str(zcase)));
        C_dualarea=compute_dualarea(C.area,C.ncell,C.nvert,C.cell_v);
        
        % Computing \zeta(u)
        CZU=zetau(CS.usol_idt,zcase);
        FZU=zetau(FS.usol_idt,zcase);

        % interpolating coarse meshes in time dimension using the solution
        % at the finest possible mesh
        ICTZ=interpolate_time(CZU,F_Ndt,meshes{imesh},meshes{nbmeshes});
        
        Cmesh=meshes{imesh};
        Fmesh=meshes{nbmeshes};
        
        parfor idt=2:F_Ndt+1
        % Interpolating coarse meshes in space dimension using the solution
        % at the finest possible mesh
            ICTSZ=interpolate_mesh(ICTZ(:,idt),FZU(:,idt),Cmesh,Fmesh);

            F_Z=FZU(:,idt);           
            errZ=F_Z-ICTSZ;
            %
            GFznorm(idt)=F_Z'*F_A*F_Z;
            Gerr_znorm(idt)=(errZ'*F_A*errZ);
            
            % 
            Fznorm(idt)=sum(F_dualarea.*abs(F_Z).^2);
            err_znorm(idt)=sum(F_dualarea.*abs(errZ).^2);
        end
        %Computing individual norms on coarse meshes 
         parfor idt=2:C_Ndt+1
             C_Z=CZU(:,idt)
             Cznorm(idt)=sum(C_dualarea.*abs(C_Z).^2);
             GCznorm(idt)=C_Z'*C_A*C_Z;
            
         end
        %For \Xi at the final time step
        CU_X=Xi(CS.usol_idt,zcase);
        FU_X=Xi(FS.usol_idt(:,F_Ndt+1),zcase);
        ICTX=interpolate_time(CU_X,F_Ndt,meshes{imesh},meshes{nbmeshes});
        ICTS_TX=interpolate_mesh(ICTX(:,F_Ndt+1),FU_X,Cmesh,Fmesh);
        %
        errX=FU_X-ICTS_TX;

        %norms in space and time for each Brownian motion
        st_eXi(bmm)=sum(F_dualarea.*abs(errX));
        Fst_Xi(bmm)=sum(F_dualarea.*FU_X);
        Cst_Xi(bmm)=sum(C_dualarea.*CU_X(:,C_Ndt+1));
        Fst_znorm(bmm)=sum(Fznorm);
        Cst_znorm(bmm)=sum(Cznorm);
        GFst_znorm(bmm)=sum(GFznorm);
        GCst_znorm(bmm)=sum(GCznorm);
        st_eznorm(bmm)=sum(err_znorm);
        Gst_eznorm(bmm)=sum(Gerr_znorm);       
    end%end Brownian motion
    %Calcuting the expectation for each mesh
    Est_eznorm(imesh)=sqrt((F_dt/nbm)*sum(st_eznorm))/sqrt((F_dt/nbm)*sum(Fst_znorm));
    GEst_eznorm(imesh)=sqrt((F_dt/nbm)*sum(Gst_eznorm))/sqrt((F_dt/nbm)*sum(GFst_znorm));
    Est_eXi(imesh)=abs((1/nbm)*sum(st_eXi))/((1/nbm)*sum(Fst_Xi));
    L2zeta(imesh)=sqrt((C_dt/nbm)*sum(Cst_znorm));
    H1zeta(imesh)=sqrt((C_dt/nbm)*sum(GCst_znorm));
    L1Xi(imesh)=(sum(Cst_Xi))/nbm;
end
    L2zeta(nbmeshes)=sqrt((F_dt/nbm)*sum(Fst_znorm));
    H1zeta(nbmeshes)=sqrt((F_dt/nbm)*sum(GFst_znorm));
    L1Xi(nbmeshes)=(sum(Fst_Xi))/nbm;

Est_eznorm
GEst_eznorm
Est_eXi
L2zeta
H1zeta
L1Xi

save(strcat('norms/T_norms_bm',num2str(nbm),'RefM',num2str(F_Ndt),'nbmeshes',num2str(nbmeshes)),'C_h','F_h','meshes','Est_eznorm','GEst_eznorm','Est_eXi','C_dof','L1Xi',"H1zeta","L2zeta");
disp(['Elapsed time is ' num2str(round(toc/60,1)) ' minutes.'])