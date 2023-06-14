
%dzeta(u)= 1 (u< 1), 0 (1\le u \le 2) and 1 (u> 2)
function dzu=dzeta(u,zcase)
dzu=zeros(size(u));
zc=zcase;
if (zc==1)
    dzu=(u<=1)+(u>2);
elseif (zc==2)
    dzu=3*u.^2;
end
end
   