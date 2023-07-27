% Compute areas of Donald dual cells

function dualarea=compute_dualarea(area,ncell,nvert,cell_v,B_indices);

dualarea=zeros(nvert,1);
for i=1:ncell
    n=3;
    multiplier=ones(3,1);
        J=find(ismember(cell_v{i}(1:3),B_indices));
        I=size(J,2);
        if I>0
            multiplier(J)=0;
            n=n-I;            
        end  
    dualarea(cell_v{i}(1:3)) = dualarea(cell_v{i}(1:3)) + (area(i)/n).*multiplier;
end



