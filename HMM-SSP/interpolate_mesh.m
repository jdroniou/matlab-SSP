function ICU=interpolate_mesh(CU,RU,Cmesh,Rmesh)

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

Cxk=gravity_centers(C.ncell,C.cell_v,C.vertex,C.area);
for i=1:size(RU,1)
    ci=C_cell(i);
    nbe=size(C.cell_e{C_cell(i)},2);
    N=zeros(nbe,2);
    N([1:nbe],:)=(C.vertex(C.cell_v{ci}([2:nbe+1]),:)-C.vertex(C.cell_v{ci}([1:nbe]),:))*[0 -1;1 0];
    % Computation of the length of the edges
    msigma=sqrt(sum(N.^2,2));
    % Normalisation of N
    N=N./[msigma msigma];
    %vector of values of edge unknowns of Ccell i
    cei=CU(C.ncell+C.cell_e{ci});
    G=sum(cei.*N.*[msigma msigma]./C.area(ci),1);
    ICU(i)=CU(ci)+G*(F_x(i,:)-Cxk(ci,:))';
      
end
