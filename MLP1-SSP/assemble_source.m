% Assemble the global diffusion matrix
function F = assemble_source(t,cell_v,ncell,nvert,area,cg,zcase,ucase)
F=zeros(nvert,1);

    for i=1:ncell
%         as1=cg(i,:)
% as2=fs(t,cg(i,:))
%  as3=F(cell_v{i}(1:3))
        F(cell_v{i}(1:3))=F(cell_v{i}(1:3))+(area(i)/3)*fs(t,cg(i,:),zcase,ucase);
    
    end
end