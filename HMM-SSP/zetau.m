%zeta(u)= u (u\le 1), 1 (1\le u \le 2) and 2 (u\ge 2)
function [Zu,DZ,D2Z]=zetau(u,zc)
Zu=zeros(size(u));
%global zc

%nargout is for the number of return values from [Zu,gradZ,D2Z]

    if (nargout>1)
        DZ=zeros(size(u));
        if (nargout>2)
            D2Z=zeros(size(u));
        end
   end
   if (zc==1)
       %Zu:=Zeta(u)={u for u<=1, 1 for 1=<u<=2, u-1 for u>=2
         Zu=u.*(u<=1)+(u>1).*(u<=2)+(u-1).*(u>2);
           if (nargout>1)
                 DZ=(u<=1)+(u>2);
               if (nargout>2)
                   %Diff Zeta''(u)=0 no-change in D2Z
               end
           end
  

   elseif (zc==2)
       Zu=u.^3;
       if (nargout>1)
           DZ=3*u.^2;
           if (nargout>2)
               D2Z=6*u;
           end
       end
   elseif (zc==3)
       Zu=u;
   end
end
