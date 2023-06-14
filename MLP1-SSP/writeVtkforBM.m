%% Writing VTK for Paraview%%%%%%%%%
% Write the averaged solution for all the Brownain motion in terms of vtk files, to be plotted by "paraview"
clear all
ucase=1;
zcase=1;
T=1;
% Select the mesh to generate the VTK files
imesh=2;
meshes={'mesh1_01.mat';'mesh1_02.mat';'mesh1_03.mat'};%'mesh1_04.mat';'mesh1_05.mat';'mesh1_06.mat'};
load(strcat('../matlab_meshes/',meshes{imesh}))
h=max(abs(diam));
Ndt=ceil(T/h^2)
%Select the number of Brownian motions
nbm=2;
Eusol_idt=zeros(nvert,Ndt+1);

parfor bmm=1:nbm
    filecontent=load(strcat('solutions/BM',num2str(bmm),'msh',meshes{imesh}(1:8),'tcuz',num2str(ucase),num2str(zcase)));
    Eusol_idt=Eusol_idt+filecontent.usol_idt;
end
Eusol_idt=(1/nbm)*Eusol_idt;
parfor idt=1:Ndt

    write_solution_vtk(Eusol_idt(:,idt),strcat('VTKout/solution',num2str(idt)),...
        ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%