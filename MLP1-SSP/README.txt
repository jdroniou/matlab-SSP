%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mass-lumped P1 code for u_t-Delta zeta(u)=f(zeta(u))*dW with Dirichlet BC
%
%    Author: Muhammad Awais Khan
%    Date: 14/06/2023
% % %
% This code is provided as a companion of the article
%   "Numerical analysis of the stochastic Stefan problem",
%
%
% Usage and modification of this code is permitted, but any scientific publication resulting from
% part of this code should cite the aforementioned article.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Run the files in the following order: 

* Step 1: Run "MAIN_MLP1_SSP.m" to save the solutions for each time step, for each Brownain motion and for each mesh.
* Step 2 [see also Note 1 below]: Run "compare_time.m" and "compare_meshes.m" to prepare the comparison between norms of solutions.
* Step 3: Run "Compute_norms.m" to compute the L2 error of: \zeta(u) and \nabla\zeta(u) and L1 error of \Xi(u).

Note 1: Step 2 only has to be run once for each mesh family; the files created in the Time-comparisons and Mesh-comparisons folders
during this step can be re-used for any solution computed on the same mesh families.
Note 2: The averaged solutions can be generated by "writeVtkforBM.m" and viewed in Paraview.
