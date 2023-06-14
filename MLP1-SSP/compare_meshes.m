clc
clear
format long;
meshes={'mesh1_01.mat';'mesh1_02.mat';'mesh1_03.mat';'mesh1_04.mat';'mesh1_05.mat';'mesh1_06.mat'};
nbmeshes=size(meshes,1);
for j=1:(nbmeshes-1);
    for k=(j+1):nbmeshes;
        % C is for coarse mesh
        % F is for finer mesh
        C=load(strcat('../matlab_meshes/',meshes{j}));
        F=load(strcat('../matlab_meshes/',meshes{k}));
        C_h=max(abs(C.diam));
        C_cell=zeros(F.nvert,1);
        Lambdas=zeros(F.nvert,3);
        C_cg=gravity_centers(C.ncell,C.cell_v,C.vertex,C.area);
        % We are looking for a triangle in the coarse mesh for the vertex s in finer
        % mesh
        parfor i=1:F.nvert
            s=F.vertex(i,:);
            % Calculating distance of s with each cell's berrycenter
            D=sqrt(sum((C_cg-s).^2,2));
            R=find(D<=2*C_h);
            for r=1:size(R,1)
                vertices=[C.vertex(C.cell_v{R(r)}(1),1) C.vertex(C.cell_v{R(r)}(2),1) C.vertex(C.cell_v{R(r)}(3),1);C.vertex(C.cell_v{R(r)}(1),2) C.vertex(C.cell_v{R(r)}(2),2) C.vertex(C.cell_v{R(r)}(3),2)];
                L=[ones(1,3);vertices]\[1;s'];

                if (all(L>=-1e-13))
                    Lambdas(i,:)=L';
                    C_cell(i)=R(r);
                end

            end
            
        end
        save(strcat('Mesh-comparisons/C',meshes{j}(1:8),'F',meshes{k}(1:8)),'C_cell','Lambdas');
% Check
%         error_location=0;
%         for i=1:F.nvert
%             error_location=max(error_location, norm(F.vertex(i,:) - Lambdas(i,1)*C.vertex(C.cell_v{C_cell(i)}(1),:) - Lambdas(i,2)*C.vertex(C.cell_v{C_cell(i)}(2),:) - Lambdas(i,3)*C.vertex(C.cell_v{C_cell(i)}(3),:) ));
%         end
%         error_location

        
    end
end
disp('Computed!')