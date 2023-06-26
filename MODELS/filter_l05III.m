% This function separates the large (x) and small (y) scales from 
% the model state (z) defined by Model III of Lorenz 2005 

function [x,y] = M_nl(z,I)

  % Parameters
  Nx    = length(z(:,1));
  alpha = (3*I^2 + 3)/(2*I^3 + 4*I);
  beta  = (2*I^2 + 1)/(I^4 + 2*I^2);
  I     = round(I);

  % Partition z into x and y
  z0 = [z;z;z];
  i = [-(I-1):I-1]';

  for m = 1:Nx
    n = Nx + m;
    x(m) = sum( (alpha - beta.*abs(i)).* z0(n+i) ) + ...
                (alpha - beta.*abs(-I)).* z0(n-I)/2  + ...
                (alpha - beta.*abs(I)).* z0(n+I)/2;
  end
  x = x';
  y = z - x;

end
