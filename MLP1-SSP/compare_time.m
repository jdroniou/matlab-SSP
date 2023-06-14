%For each fine time t, find the coarse time T_n such that T_n<=t<T_{n+1} (it's simply n=ceil(t/dt),
% if dt is the coarse time step, I think), then write t=aT_n+(1-a)T_{n+1}
% (so a=(T_{n+1}-t)/(T_{n+1}-T_n)).
% Compute the coarse solution at time t by taking "a u(T_n)+(1-a)u(T_{n+1})".

clc
clear
format long;
ucase=1;
zcase=1;
T=1;
%% Select meshes
meshes={'mesh1_01.mat';'mesh1_02.mat';'mesh1_03.mat';'mesh1_04.mat';'mesh1_05.mat';'mesh1_06.mat'};
nbmeshes=size(meshes,1);
for j=1:(nbmeshes-1);
    for k=(j+1):nbmeshes;
        clear n a
        C=load(strcat('../matlab_meshes/',meshes{j}));
        F=load(strcat('../matlab_meshes/',meshes{k}));
        C_h=max(abs(C.diam));
        F_h=max(abs(F.diam));
        % Number of time steps are computed by assuming dt=h^2
        C_Ndt=ceil(T/C_h^2);
        F_Ndt=ceil(T/F_h^2);
        C_dt=T/C_Ndt;
        F_dt=T/F_Ndt;
        n=zeros(F_Ndt,1);
        V_T=zeros(F_Ndt,1);
        a=zeros(F_Ndt,1);
        Ft=zeros(F_Ndt,1);
        V_t=F_dt*[1:F_Ndt]';
        n=ceil(V_t/C_dt);
        T_n=(n-1)*C_dt;
        T_n1=n*C_dt;
        a=(T_n1-V_t)./(T_n1-T_n);

         save(strcat('Time-comparisons/C',meshes{j}(1:8),'F',meshes{k}(1:8)),'n','a');
       
    end
end
disp('Computed!')


