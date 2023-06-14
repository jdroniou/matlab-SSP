% Write a .vtk file, readable by paraview, containing the
% mesh and the solution 
%
% u should be a vector of size ncell+nedge, with the cell values first
%
function out=write_solution_vtk_P1(u,namefile,ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex)

%
% We write the vtk file
%

completename=strcat(namefile,'.vtk');
fid=fopen(completename,'w');


fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'Grille\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET POLYDATA\n');
fprintf(fid,'POINTS %d float\n',nvert);
for i=1:nvert
	fprintf(fid,'%E %E %E\n',vertex(i,1),vertex(i,2),u(i));
end
sumvert=0;
for i=1:ncell
	sumvert=sumvert+size(cell_v{i},2)-1;
end
fprintf(fid,'POLYGONS %d %d\n',ncell,sumvert+ncell);
for i=1:ncell
% Should work for all type of mesh
		fprintf(fid,'%d\n',size(cell_v{i},2)-1,cell_v{i}(1:size(cell_v{i},2)-1)-1);
end
fprintf(fid,'CELL_DATA %d\n',ncell);
fprintf(fid,'SCALARS sol float 1\n');
fprintf(fid,'LOOKUP_TABLE default\n');
%Compute of the gravity center of the mesh
u_g=zeros(ncell,1);
for i=1:ncell
    u_g(i)=(u(cell_v{i}(1))+u(cell_v{i}(2))+u(cell_v{i}(3)))/3;
end

for i=1:ncell
	fprintf(fid,'%E\n',u_g(i));
end

fclose(fid);




