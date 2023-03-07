% This function integrates the Lorenz-05 equations forward in time 
% using the matlab 'ode45' solver.
%
%  INPUT 
%         x: initial state (Nx x 1)
%        dt: time step
%         T: integration time
%         F: forcing term
%
%  OUTPUT 
%         x: model state integrated from initial state (Nx x 1)

function z = M_nl(z,dt,T,K,I,b,c,F)

  % Parameters
  Nx    = length(z(:,1));
  K     = round(K);
  I     = round(I);
  alpha = (3*I^2 + 3)/(2*I^3 + 4*I);
  beta  = (2*I^2 + 1)/(I^4 + 2*I^2);

  % Time integration
  tspan=[0,T]*dt;
  options=' ';
%  options = odeset('RelTol',1e-5);
  [~,z] = ode45(@dzt, tspan, z(:,1)',options,K,I,b,c,F,alpha,beta,Nx);
  z = z(end,:)';

  % Function for dz
  function [z1] = dzt(~,z0,K,I,b,c,F,alpha,beta,Nx)

    % Partition z into x and y
    z0 = [z0;z0;z0];
    i = [-(I-1):I-1]';

    if I == 1

      x0 = z0;

    else

      for m = 1:Nx
        n = Nx + m;
        x0(m,1) = sum( (alpha - beta.*abs(i)).* z0(n+i) ) + ...
                   (alpha - beta.*abs(-I)).* z0(n-I)/2  + ...
                   (alpha - beta.*abs(I)).* z0(n+I)/2;
      end

      y0 = z0(Nx+1:2*Nx) - x0;

    end

    % Create buffer zones for x0 and y0
    x0 = [x0;x0;x0]; 
    if I > 1, y0 = [y0;y0;y0]; end
    w  = zeros(size(x0));

    % Indices for spatial smoothing part
    J = floor(K/2);
    j = [-(J-1):J-1]';

    if mod(K,2) == 0
      norm = 1/2;
    else
      norm = 1;
    end

    % Weight vectors for XX calculations;
    % see (9) -- (10) of Lorenz (2005)
    for m = Nx-2*K:2*Nx+2*K
      w(m) = ( sum( x0(m-j) ) + (x0(m-J) + x0(m+J))*norm ) / K;
    end

    % XX term
    for m = 1:Nx
      n = Nx + m;
      xx(m,1) = -w(n-2*K)*w(n-K) + ( sum( w(n-K+j).*x0(n+K+j) )+...
         ( w(n-K-J).*x0(n+K-J) + w(n-K+J).*x0(n+K+J) )*norm ) / K;
    end

    % Indices for selecting states on period grid
    i1 = Nx+[-1:Nx-2]; i2 = Nx+[0:Nx-1];
    i3 = Nx+[1:Nx]; i4 = Nx+[2:Nx+1];

    if I > 1

      % YY and YX terms
      yy = -y0(i1,1).*y0(i2,1) + y0(i2,1).*y0(i4,1);
      yx = -y0(i1,1).*x0(i2,1) + y0(i2,1).*x0(i4,1);
    
      % Final advance step
      z1 = xx + b^2*yy + c*yx - x0(i3,1) - b*y0(i3,1) + F;

    else

      % Final advance step
      z1 = xx - x0(i3,1) + F;

    end

  end

end
