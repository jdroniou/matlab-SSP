% Assemble the global diffusion matrix
function F = assemble_source(zc,idt,cell_v,ncell,nvert,area,dW,Xprev)
F=zeros(nvert,1);
zcase=zc;
    for i=1:ncell
        F(cell_v{i}(1:3))=F(cell_v{i}(1:3))+(area(i)/3)*dW(idt)*sqrt(Xi(Xprev(cell_v{i}(1:3)),zcase));
    
    end
end