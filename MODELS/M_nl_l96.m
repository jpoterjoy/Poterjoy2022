% This function integrates the Lorenz-96 equations forward in time 
% using the fourth order Runge-Kutta scheme
%
%  INPUT 
%         x: initial state (Nx x 1)
%        dt: time step
%         T: integration time
%         F: forcing term
%
%  OUTPUT 
%         x: model state integrated from initial state (Nx x 1)

function x = M_nl(x,dt,T,F)

  Nx = length(x(:,1));

  % Create function for dx/dt
  function x1 = dxt(x0,F)

    % Create buffer zones on x for periodic domain
    x0 = [x0(Nx-1,1); x0(Nx,1); x0(:,1); x0(1,1)];

    % Place variables in vectors
    y1 = x0(1:end-3,1); y2 = x0(2:end-2,1);
    y3 = x0(3:end-1,1); y4 = x0(4:end,1);   
    x1 = -y2.*(y1 - y4) - y3 + F;

  end

  for t = 1:T

    % Find RK coefficients
    k1 = dt*dxt( x       , F );
    k2 = dt*dxt( x + k1/2, F );
    k3 = dt*dxt( x + k2/2, F );
    k4 = dt*dxt( x + k3  , F );

    % Update variables
    x(:,1) = x(:,1) + k1/6 + k2/3 + k3/3 + k4/6;

  end

end
