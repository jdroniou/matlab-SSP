%integral of zeta(u)= u (u\le 1), 1 (1\le u \le 2) and 2 (u\ge 2)
%xi = int(piecewise(u<=1,u,u>1 & u <2,1,u >2,u-1),0,upperlimit);

function xi=Xi(u,zc)%Xi(u)=int(zeta,0,u)
zcase=zc;
if (zcase==1)
    xi=(0.5*u.^2).*(u<=1) + u.*(u>1).*(u<=2) + 0.5*u.*(u-2).*(u>2);
elseif (zcase==2)
    xi=(0.25*(u.^4))';
end
end