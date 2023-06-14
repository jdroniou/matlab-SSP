% Transform between the real function u, and the vector of unknowns X used in computations
	%		m>1: X is the solution 'u' in the cells (1:ncell) and '|u|^(m-1)u' on the edges (ncell+1:ncell+nedge)
  %   m<1: X is |u|^(m-1)u everywhere

function a=transform(b,ncell,direction)

global m;
sizeb = size(b,1);
a=zeros(size(b));
sb = sign(b);

if (direction=='to_unknowns')
  % b is the real function, a are the unknowns
   a(1:ncell) = b(1:ncell);
   a(ncell+1:sizeb) = zetau(b(ncell+1:sizeb));%abs(b(ncell+1:sizeb)).^m.*sb(ncell+1:sizeb);

end


