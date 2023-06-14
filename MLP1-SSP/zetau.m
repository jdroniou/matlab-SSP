function [Zu,DZ,D2Z]=zetau(u,zcase)
zc=zcase;
sz=size(u);
Zu=zeros(sz);
%nargout is for the number of return values from [Zu,gradZ,D2Z]

if (nargout>1)
    DZ=zeros(sz);
    if (nargout>2)
        D2Z=zeros(sz);
    end
end
if (zc==1)
    %Zu:=Zeta(u)={u for u<=1, 1 for 1=<u<=2, u-1 for u>=2
    for i=1:sz
        Zu=u.*(u<=1)+(u>1).*(u<=2)+(u-1).*(u>2);
        if (nargout>1)
            DZ=(u<=1)+(u>2);
            if (nargout>2)
                %Diff Zeta''(u)=0 no-change in D2Z
            end
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
end

end