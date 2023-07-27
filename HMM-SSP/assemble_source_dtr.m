% Computing source term
function F = assemble_source(t,Area_CE,all_cg,zc,cell_e,ncell,uc)
n=size(Area_CE,1);
F=zeros(n,1);
for i=1:n
    
    F(i)=Area_CE(i)*fs(t,all_cg(i,:),zc,uc);
    % nbe=size(cell_e{i},2);
    % for j=1:nbe
    %     k=cell_e{i}(j);
    %     F(k)=F(k)+Area_CE(k)*fs(t,all_cg(k,:),zc);
    % end
end
end
