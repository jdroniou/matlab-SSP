% Computing source term
function F = assemble_source(zc,idt,area_cells,area_edges,dW,Xprev)
zcase=zc;
F=dW(idt)*[area_cells; area_edges].*sqrt(Xi(Xprev,zcase));
end