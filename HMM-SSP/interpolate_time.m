function ICT=interpolate_time(CU,F_Ndt,Cmesh,Rmesh)

tc=load(strcat('Time-comparisons/C',Cmesh(1:8),'F',Rmesh(1:8)));
a=tc.a;
n=tc.n;
ICT=zeros(size(CU,1),F_Ndt+1);
% size(ICT)
ICT(:,1)=CU(:,1);
for idt=2:F_Ndt+1
  
        ICT(:,idt)=a(n(idt-1))*CU(:,n(idt-1))+(1-a(n(idt-1)))*CU(:,n(idt-1)+1);
  
end
