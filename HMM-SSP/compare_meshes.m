clc
clear
format longE;
%Hexagonal meshes
%meshes={'hexa1_01.mat';'hexa1_02.mat';'hexa1_03.mat';'hexa1_04.mat';'hexa1_05.mat'};
%Triangular meshes
 meshes={'mesh1_01.mat';'mesh1_02.mat';'mesh1_03.mat';'mesh1_04.mat';'mesh1_05.mat'};
nbmeshes=size(meshes,1);
for j=1:(nbmeshes-1);
    for k=(j+1):nbmeshes;
        fprintf('%s and %s\n', meshes{j},meshes{k})
        % C is for coarse mesh
        % F is for finer mesh
        C=load(strcat('../matlab_meshes/',meshes{j}));
        F=load(strcat('../matlab_meshes/',meshes{k}));
        C_h=max(abs(C.diam));
        C_cell=zeros(F.ncell+F.nedge,1);
        Lambdas=zeros(F.nvert,3);
        F_cg=gravity_centers(F.ncell,F.cell_v,F.vertex,F.area);
        C_cg=gravity_centers(C.ncell,C.cell_v,C.vertex,C.area);
        % Vector of vertices of unknowns of finer mesh
        F_un=zeros(F.ncell+F.nedge,2);
        for i=1:F.ncell
            nbe = size(F.cell_e{i},2);
            F_un(i,:)=F_cg(i,:);
            F_un(F.ncell+F.cell_e{i},:)=(F.vertex(F.cell_v{i}(1:nbe),:)+F.vertex(F.cell_v{i}(2:nbe+1),:))/2;
        end
  
        % We are looking for an element in the coarse mesh for the vertex s in the finer
        % mesh
        for i=1:(F.ncell+F.nedge)
            s=F_un(i,:);
            % Calculating the distance of s with each cell's berrycenter
            D=sqrt(sum((C_cg-s).^2,2));
            R=find(D<=(1+1e-3)*C_h);
%             pos=0;
            for r=1:size(R,1)
                nbe = size(C.cell_e{R(r)},2);
                % Center of edges of coarse cell i
                xec=(C.vertex(C.cell_v{R(r)}(1:nbe),:)+C.vertex(C.cell_v{R(r)}(2:nbe+1),:))/2;
                %
                N=zeros(nbe,2);
                N([1:nbe],:)=(C.vertex(C.cell_v{R(r)}(2:nbe+1),:)-C.vertex(C.cell_v{R(r)}(1:nbe),:))*[0 -1;1 0];
                % Computation of the length of the edges
                msigma=sqrt(sum(N.^2,2));
                dv=s-xec;
                %mxxec=sqrt(sum(dv.^2,2));
              

                % Normalisation
                N=N./[msigma msigma];
                %dv=dv./[mxxec mxxec];
                dp=sum((dv).*N,2);

                if (all(dp<=1e-13))

                    C_cell(i)=R(r);

                end

            end

        end
        % Testing that if any enrty in the indexed vector C_cell is zero
       if (size(find(~C_cell),1)==0)
          save(strcat('Mesh-comparisons/C',meshes{j}(1:8),'F',meshes{k}(1:8)),'C_cell');
       else
           size(find(~C_cell))
           fprintf('Error at %s and %s\n',meshes{j},meshes{k})
           find(~C_cell)
           return
       end


    end
end
disp('Computed!');