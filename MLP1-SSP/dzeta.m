
%dzeta(u)= 1 (u< 1), 0 (1\le u \le 2) and 1 (u> 2)
function dzu=dzeta(u,zcase)
dzu=zeros(size(u));
% global zcase
if (zcase==1)
    for i=1:size(u)
        dzu=(u<=1)+(u>2);
    end
elseif (zcase==2)
    dzu=3*u.^2;
elseif (zcase==3)
    dzu=1;
end
end
