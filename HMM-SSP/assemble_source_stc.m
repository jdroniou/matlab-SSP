% Computing source term
function F = assemble_source_stc(zc,idt,Area_CE,dW,Xprev)

F=dW(idt)*Area_CE.*sqrt(Xi(Xprev,zc));
end