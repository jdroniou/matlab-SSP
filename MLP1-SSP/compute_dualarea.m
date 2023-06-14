% Compute areas of Donald dual cells

function dualarea=compute_dualarea(area,ncell,nvert,cell_v);

dualarea=zeros(nvert,1);
for i=1:ncell
  for j=1:3
    dualarea(cell_v{i}(j)) = dualarea(cell_v{i}(j)) + area(i)/3;
	end
end



