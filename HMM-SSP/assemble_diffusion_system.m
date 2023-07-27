% Assemble the global diffusion matrix
function [DiffMat,source] = assemble_diffusion_system(cell_v,cell_n,cell_e,ncell,nedge,vertex,area,center,cg);


	source=zeros(ncell+nedge,1);

	%% Initialise vectors for sparse matrix
	% Evaluate number of non-zeros entries (cf below how many times we do "pos=pos+1")
	nz=0;
	for i=1:ncell
		nedge_i=size(cell_e{i},2);
		nz = nz + nedge_i * (1+4*nedge_i);
	end
	IA=zeros(nz,1);
	JA=zeros(nz,1);
	VA=zeros(nz,1);

	% "pos"=position inside the vectors IA, JA, VA that store the entries of A
	pos=0;
	% We fill in the matrix as per the lecture notes of the AMSI summer school 2016
	%
	%% Loop over cells
	for i=1:ncell
		% Compute local matrices
		W=local_flux_matrix(vertex(cell_v{i},:),area(i),center(i,:),cg(i,:));

		% Loop over edges of cell
		for jj=1:size(cell_e{i},2)
			j=cell_e{i}(jj);
% 			% Dirichlet boundary condition
% 			if (cell_n{i}(jj)==0)
% 				pos=pos+1;
% 				IA(pos)=ncell+j;
% 				JA(pos)=ncell+j;
% 				VA(pos)=1;
% 			end
			
			% Inner loop over edges of cell
			for kk=1:size(cell_e{i},2)
				k=cell_e{i}(kk);

				% Entry cell-cell
				pos=pos+1;
				IA(pos)=i;
				JA(pos)=i;
				VA(pos)=W(jj,kk);

				% Entry cell-edge k
				pos=pos+1;
				IA(pos)=i;
				JA(pos)=ncell+k;
				VA(pos)=-W(jj,kk);

				% If j is not a Dirichlet boundary edge
% 				if (cell_n{i}(jj)~=0)
					% Entry edge j-cell
					pos=pos+1;
					IA(pos)=ncell+j;
					JA(pos)=i;
					VA(pos)=-W(jj,kk);

					% Entry edge j-edge k
					pos=pos+1;
					IA(pos)=ncell+j;
					JA(pos)=ncell+k;
					VA(pos)=W(jj,kk);
% 				end
			end;
		end;	
		

	end;

	%% Creation of the sparse matrix
	DiffMat=sparse(IA(1:pos),JA(1:pos),VA(1:pos),ncell+nedge,ncell+nedge);

