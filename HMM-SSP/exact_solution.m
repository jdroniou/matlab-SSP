% Exact solution: only used to compute intial and boundary conditions.
%This solution is taken from the article "Uniform-in-time convergence of numerical
%schemes for Richards% and Stefan"s models" by
%Jérôme Droniou, Robert Eymard and Cindy Guichard.
function [u,gradu,D2u]=exact_solution(t,z)

global uc;

x=z(:,1);
y=z(:,2);
%u=zeros(size(z,1),1);
%nargout is for the number of return values from [u,ut,gradu,D2u]
if (nargout>1)
    gradu=zeros(size(z,1),2);
    if (nargout>2)
        D2u=zeros(2,2,size(z,1));
    end
end

if (uc==1)
    u=(2*exp(t-x)).*(x<=t)+(exp(t-x)).*(x>t);
    if (nargout>1)
        %compute div(grad zeta)
        gradu(:,1)=-u;
        %gradz(:,2)=0;
        % No change to the 2nd component of the gradient, remains 0
        if (nargout>2)
            % D2z, all others remain 0
            D2u(1,1,:)=u;
            % 			    D2u(1,2,:)=0;
            % 			    D2u(2,1,:)=0;
            % 			    D2u(2,2,:)=0;
        end
    end
  

elseif (uc==2)
    % u(x,y)=cos(t)sin(pi x)sin(pi y)
    u=1+(cos(t)*(sin(pi*x).*sin(pi*y)))';
%     if (nargout>1)
%         ut(:,1)=-sin(t)*(sin(pi*x).*sin(pi*y));
%         if (nargout>2)
%             gradu(:,1)=cos(t)*(pi*cos(pi*x).*sin(pi*y));
%             gradu(:,2)=cos(t)*(pi*cos(pi*y).*sin(pi*x));
%             if (nargout>3)
%                 D2u(1,1,:)=cos(t)*(-pi^2*sin(pi*x).*sin(pi*y));
%                 D2u(1,2,:)=cos(t)*(pi^2*cos(pi*x).*cos(pi*y));
%                 D2u(2,1,:)=D2u(1,2,:);
%                 D2u(2,2,:)=D2u(1,1,:);
%             end
%         end
%     end

end
end
