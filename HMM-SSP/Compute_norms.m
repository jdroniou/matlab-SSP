%In the following code, we will be calculating error in relavant norms for
% the reconstructed functions of \zeta(u), \nabla\zeta(u) and \Xi(u)
clear
format long;
tic
ucase=1;
zcase=1;
T=1;
nbm=10;
%r is a parameter to disturbute a portion of Mass from the cells center to the edges
r=0.5;
%% Select meshes
%Hexagonal meshes
%meshes={'hexa1_01.mat';'hexa1_02.mat';'hexa1_03.mat'};%'hexa1_04.mat';'hexa1_05.mat'};
%Triangular meshes
meshes={'mesh1_01.mat';'mesh1_02.mat';'mesh1_03.mat';'mesh1_04.mat'};%'mesh1_05.mat'};
nbmeshes=size(meshes,1);
% F is for finer mesh
F=load(strcat('../matlab_meshes/',meshes{nbmeshes}));
Fcg=gravity_centers(F.ncell,F.cell_v,F.vertex,F.area);
F_A=assemble_diffusion_system(F.cell_v,F.cell_n,F.cell_e,F.ncell,F.nedge,F.vertex,F.area,F.center,Fcg);
F_h=max(abs(F.diam));
F_Ndt=2*round(0.5*T/F_h^2);
F_dt=T/F_Ndt;
% Separating boundary and interior unknowns
 F_aec=zeros(F.nedge,2);
    F_fbe=zeros(F.nedge,1);
for i=1:F.ncell
%         Fnbe=size(F.cell_e{i},2);
%         aec(F.cell_e{i},:)=(F.vertex(F.cell_v{i}(1:Fnbe),:)+F.vertex(F.cell_v{i}(2:Fnbe+1),:))/2;
        I=find(F.cell_n{i}==0);
        if (size(I,2)>0)
            % Fixing it in order to find the boundary edges
            F_fbe(F.cell_e{i}(I))=1;
        end;
    end

    F_BE_indices=find(F_fbe);
    F_IE_indices=find(~F_fbe);
    F_int_all=[1:F.ncell F.ncell+F_IE_indices']';
    % Global boundary indices
    F_GBE_indices=F.ncell+F_BE_indices;
%% Distributing mass among edges
Farea_edges=zeros(F.nedge,1);
for i=1:F.ncell
    Fnbe=size(F.cell_e{i},2);
    Frs=ones(Fnbe,1);
     rs=ones(Fnbe,1);
        FJ=find(ismember(F.cell_e{i},F_BE_indices));
        FI=size(FJ,2);
        if FI>0
            Frs(FJ)=0;
            Fnbe=Fnbe-FI;
        end   
    Farea_edges(F.cell_e{i})=Farea_edges(F.cell_e{i})+((1-r)*F.area(i)).*Frs/Fnbe;
end;
Farea_cells=r.*F.area;
F_Area_CE=[Farea_cells;Farea_edges];


Fst_znorm=zeros(nbm,1);
GFst_znorm=zeros(nbm,1);
Fst_Xi=zeros(nbm,1);
Fznorm=zeros(F_Ndt+1,1);

fid = fopen('results.txt','w');

Est_eznorm=zeros(nbmeshes-1,1);
GEst_eznorm=zeros(nbmeshes-1,1);
Est_eXi=zeros(nbmeshes-1,1);
L2zeta=zeros(nbmeshes-1,1);
H1zeta=zeros(nbmeshes-1,1);
L1Xi=zeros(nbmeshes-1,1);
C_dof=zeros(nbmeshes-1,1);
C_h=zeros(nbmeshes-1,1);

for imesh=1:nbmeshes-1
    % C is for coarse mesh
    C=load(strcat('../matlab_meshes/',meshes{imesh}));
    str = sprintf('Loading %s and %s\n',meshes{imesh},meshes{nbmeshes});
    forkprint(fid,str);
    Ccg=gravity_centers(C.ncell,C.cell_v,C.vertex,C.area);
    C_A=assemble_diffusion_system(C.cell_v,C.cell_n,C.cell_e,C.ncell,C.nedge,C.vertex,C.area,C.center,Ccg);
      % Separating boundary and interior unknowns
     C_aec=zeros(C.nedge,2);
    C_fbe=zeros(C.nedge,1);
for i=1:C.ncell
%         Cnbe=size(C.cell_e{i},2);
%         aec(C.cell_e{i},:)=(C.vertex(C.cell_v{i}(1:Cnbe),:)+C.vertex(C.cell_v{i}(2:Cnbe+1),:))/2;
        I=find(C.cell_n{i}==0);
        if (size(I,2)>0)
            % Fixing it in order to find the boundary edges
            C_fbe(C.cell_e{i}(I))=1;
        end;
    end

    C_BE_indices=find(C_fbe);
    C_IE_indices=find(~C_fbe);
    C_int_all=[1:C.ncell C.ncell+C_IE_indices']';
    % Global boundary indices
    C_GBE_indices=C.ncell+C_BE_indices;
    Carea_edges=zeros(C.nedge,1);
    for i=1:C.ncell
        Cnbe=size(C.cell_e{i},2);
        Crs=ones(Cnbe,1);
        CJ=find(ismember(C.cell_e{i},C_BE_indices));
        CI=size(CJ,2);
        if CI>0
            Crs(CJ)=0;
            Cnbe=Cnbe-CI;
        end
        Carea_edges(C.cell_e{i})=Carea_edges(C.cell_e{i})+((1-r)*C.area(i)).*Crs/Cnbe;
    end;
    Carea_cells=r.*C.area;
    C_Area_CE=[Carea_cells;Carea_edges];
    C_h(imesh)=max(abs(C.diam));
    %for imesh=1:nbmeshes
    C_Ndt=2*round(0.5*T/C_h(imesh)^2);
    C_dt=T/C_Ndt;
    C_dof(imesh)=C.nedge;
  

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

        CZU=zetau(CS.usol_idt,zcase);
        FZU=zetau(FS.usol_idt,zcase);

        Csol=CS.usol_idt;
        Fsol=FS.usol_idt;
        CZU=Csol;
        FZU=Fsol;
        CZU(C_int_all,:)=zetau(Csol(C_int_all,:),zcase);
        FZU(F_int_all,:)=zetau(Fsol(F_int_all,:),zcase);


        %% Time interpolation
        % interpolating coarse meshes in time dimension using the solution
        % at the finest possible mesh
        ICTZ=interpolate_time(CZU,F_Ndt,meshes{imesh},meshes{nbmeshes});
        %% Space interpolation
        Cmesh=meshes{imesh};
        Fmesh=meshes{nbmeshes};
        parfor idt=2:F_Ndt+1
            % Interpolating coarse meshes in space dimension using the solution
            % at the finest possible mesh
            ICTSZ=interpolate_mesh(ICTZ(:,idt),FZU(:,idt),Cmesh,Fmesh);

            % Computing for L2 norm of \zeta(u_R)-\zeta(u_m) in Thetha, (0,T) and Omega
            F_Z=FZU(:,idt);
            errZ=(FZU(:,idt)-ICTSZ);
            % Computing L2 norm of Grad-zeta
            GFznorm(idt)=F_Z'*F_A*F_Z;
            Gerr_znorm(idt)=(errZ'*F_A*errZ);
            % Computing L2 norm of zeta
            Fznorm(idt)=sum(F_Area_CE.*abs(F_Z).^2);
            err_znorm(idt)=sum(F_Area_CE.*abs(errZ).^2);
        end
        %Computing individual norms on coarse meshes
        parfor idt=2:C_Ndt+1
            C_Z=CZU(:,idt)
            Cznorm(idt)=sum(C_Area_CE.*abs(C_Z).^2);
            GCznorm(idt)=C_Z'*C_A*C_Z;

        end

        CU_X=Xi(CS.usol_idt,zcase);
        FUX=Xi(FS.usol_idt,zcase);
        ICTX=interpolate_time(CU_X,F_Ndt,meshes{imesh},meshes{nbmeshes});
        ICTS_TX=interpolate_meshXi(ICTX(:,F_Ndt+1),FUX(:,F_Ndt+1),Cmesh,Fmesh);
        errX=FUX(:,F_Ndt+1)-ICTS_TX;
        
%         st_eXi(bmm)=(sum(F_Area_CE(F_int_all).*abs(errX(F_int_all))));
%         Fst_Xi(bmm)=sum(F_Area_CE(F_int_all).*FUX(F_int_all,F_Ndt+1));
%         Cst_Xi(bmm)=sum(C_Area_CE(C_int_all).*CU_X(C_int_all,C_Ndt+1));

        st_eXi(bmm)=(sum(F_Area_CE.*abs(errX)));
        Fst_Xi(bmm)=sum(F_Area_CE.*FUX(:,F_Ndt+1));
        Cst_Xi(bmm)=sum(C_Area_CE.*CU_X(:,C_Ndt+1));


        Fst_znorm(bmm)=sum(Fznorm);
        Cst_znorm(bmm)=sum(Cznorm);
        GFst_znorm(bmm)=sum(GFznorm);
        GCst_znorm(bmm)=sum(GCznorm);
        st_eznorm(bmm)=sum(err_znorm);
        Gst_eznorm(bmm)=sum(Gerr_znorm);

    end%end Brownian motion
    Est_eznorm(imesh)=sqrt((F_dt/nbm)*sum(st_eznorm))/sqrt((F_dt/nbm)*sum(Fst_znorm));
    GEst_eznorm(imesh)=sqrt((F_dt/nbm)*sum(Gst_eznorm))/sqrt((F_dt/nbm)*sum(GFst_znorm));
    Est_eXi(imesh)=abs(sum(st_eXi))/(sum(Fst_Xi));
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


save(strcat('norms/norms_bm',num2str(nbm),meshes{nbmeshes}(1:4),'RefM',num2str(F_Ndt),'nbmeshes',num2str(nbmeshes)),'C_h','F_h','meshes','Est_eznorm','GEst_eznorm','Est_eXi','C_dof','L1Xi',"H1zeta","L2zeta");
disp(['Elapsed time is ' num2str(round(toc/60,1)) ' minutes.'])