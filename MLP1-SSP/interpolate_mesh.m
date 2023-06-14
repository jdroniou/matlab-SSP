function ICU=interpolate_mesh(CU,RU,Cmesh,Rmesh)
C=load(strcat('../matlab_meshes/',Cmesh));
%R=load(strcat('../matlab_meshes/',Rmesh));
mc=load(strcat('Mesh-comparisons/C',Cmesh(1:8),'F',Rmesh(1:8)));
C_cell=mc.C_cell;
ICU=zeros(size(RU,1),1);
Lambdas=mc.Lambdas;
for i=1:size(RU,1)
    ICU(i)=Lambdas(i,1)*CU(C.cell_v{C_cell(i)}(1)) + Lambdas(i,2)*CU(C.cell_v{C_cell(i)}(2)) + Lambdas(i,3)*CU(C.cell_v{C_cell(i)}(3));
end
