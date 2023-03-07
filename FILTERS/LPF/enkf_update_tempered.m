% Function for performing enkf update using Whitaker and Hamill filter with
% Anderson adaptive state-space inflation
%
%  INPUT 
%           x: prior ensemble (Nx x Ne)
%          hx: obs-space prior ensemble (Ny x Ne)
%           y: observation vector  (Ny x 1)
%       var_y: obs error variance
%          HC: matrix determining localization between obs- and model-space  (Ny x Nx)
%         HCH: matrix determining localization between obs- and obs-space  (Ny x Ny)
%    inf_flag: flag for inflation option
%         inf: prior mean of inf coefficients (Nx x 1)
%     var_inf: prior variance of inf coefficients
%       gamma: relaxation parameter for RTPS 
%
%  OUTPUT 
%          xm: posterior mean (Nx x 1)
%           x: posterior ensemble (Nx x Ne)
%         inf: posterior mean of inf coefficients (Nx x 1)
%      e_flag: error flag

function [xm,x,e_flag] = enkf_update(x,hx,y,var_y,HC,HCH,gamma,obs_infl)

% Get array dimensions
[Nx,Ne] = size(x);
Ny = length(y);

% Mean and perturbation states
xm = mean(x')';
xp = x - xm;
hxm = mean(hx')';
hxp = hx - hxm;
xpo = xp;

for i = 1:Ny

  if obs_infl(i) > 99999, continue, end

  % Innovations
  d = y(i) - hxm(i);
  hxo = hxp(i,:);
  var_den = hxo*hxo'/(Ne-1) + var_y*obs_infl(i);

  % --- State-space update ---

  % Calculate localized gain matrix
  P = xp*hxo'/(Ne-1);
  P = P.*HC(i,:)';
  K = P/var_den;

  % Update mean
  xm = xm + K*d';

  % Update perturbations 
  beta = 1/(1 + sqrt( var_y*obs_infl(i)/var_den) );
  xp = xp - beta*K*hxo;

  % Error check
  if sum(sum(isnan(xm))) > 0
    e_flag = 1; 
    x = xp;
    return
  else
    e_flag = 0;
  end

  % --- Obs-space update ---

  % Calculate localized gain matrix
  P = hxp*hxo'/(Ne-1);
  P = P.*HCH(i,:)';
  K = P/var_den;

  % Update mean
  hxm = hxm + K*d';

  % Update perturbations 
  beta = 1/(1 + sqrt( var_y*obs_infl(i)/var_den) );
  hxp = hxp - beta*K*hxo;

end

% Apply RTPS
%  xp = xp*(1-gamma) + xpo*gamma;
v1 = sqrt(var(xpo')');
v2 = sqrt(var(xp')');
xp = xp.*( gamma*(v1 - v2)./v2 + 1 );

% Add perturbations to mean
x = xm + xp;

if sum(isnan(xm)) > 0
  xm = nan;
  x = nan
end
