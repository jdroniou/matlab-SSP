function ICU=interpolate_meshXi(CU,RU,Cmesh,Rmesh)

C=load(strcat('../matlab_meshes/',Cmesh));
F=load(strcat('../matlab_meshes/',Rmesh));
Fcg=gravity_centers(F.ncell,F.cell_v,F.vertex,F.area);
 F_x=zeros(F.ncell+F.nedge,2);
        for i=1:F.ncell
            nbe = size(F.cell_e{i},2);
            F_x(i,:)=Fcg(i,:);
            F_x(F.ncell+F.cell_e{i},:)=(F.vertex(F.cell_v{i}(1:nbe),:)+F.vertex(F.cell_v{i}(2:nbe+1),:))/2;
        end
mc=load(strcat('Mesh-comparisons/C',Cmesh(1:8),'F',Rmesh(1:8)));
C_cell=mc.C_cell;
ICU=zeros(size(RU,1),1);

%Cxk=gravity_centers(C.ncell,C.cell_v,C.vertex,C.area);
for i=1:size(RU,1)
    ci=C_cell(i);
    ICU(i)=CU(ci);
      
end
