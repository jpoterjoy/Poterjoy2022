% Model integration using Lorenz (1996) 3-variable model
%
% function x = M_nl(x,dt,T,s,r,b)

function x = M_nl(x,dt,T,s,r,b)

x2 = x;
% Forward integration
for t = 1:T
  x(1) = x2(1) + dt*s*(x2(2)-x2(1));
  x(2) = x2(2) + dt * ( r*x2(1) - x2(2) - x2(1)*x2(3));
  x(3) = x2(3) + dt * ( x2(1)*x2(2) - b*x2(3) );
  x2 = x;
end
