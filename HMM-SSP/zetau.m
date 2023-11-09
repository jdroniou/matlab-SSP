%zeta(u)= u (u\le 1), 1 (1\le u \le 2) and 2 (u\ge 2)
function [Zu,DZ,D2Z]=zetau(u,zcase)
Zu=zeros(size(u));
% global zcase

%nargout is for the number of return values from [Zu,gradZ,D2Z]

    if (nargout>1)
        DZ=zeros(size(u));
        if (nargout>2)
            D2Z=zeros(size(u));
        end
   end
      if (zcase==1)
       %Zu:=Zeta(u)={u for u<=1, 1 for 1=<u<=2, u-1 for u>=2
         Zu=u.*(u<=1)+(u>1).*(u<=2)+(u-1).*(u>2);
           if (nargout>1)
                 DZ=(u<1)+(u>2);
               if (nargout>2)
                   %Diff Zeta''(u)=0 no-change in D2Z
               end
           end
   elseif (zcase==2)
       %Zu:=Zeta(u)={u for u<=0, 0 for 0=<u<=1, u-1 for u>=1
         Zu=u.*(u<=0)+(0).*(u>=0).*(u<=1)+(u-1).*(u>=1);
           if (nargout>1)
                 DZ=(u<0)+(u>1);
               if (nargout>2)
                   %Diff Zeta''(u)=0 no-change in D2Z
               end
           end

   elseif (zcase==3)
       Zu=u.^2;
       if (nargout>1)
           DZ=2*u;
           if (nargout>2)
               D2Z=2;
           end
       end
   elseif (zcase==4)
       Zu=u;
        if (nargout>1)
           DZ=1;
           if (nargout>2)
               D2Z=0;
           end
       end
   end
end
