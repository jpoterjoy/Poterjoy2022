% Function for performing the particle filter update step
%
% * This version performs tempering steps based on regularization coefficients
%
%  INPUT 
%           x: prior particles (Nx x Ne)
%          hx: obs-space prior particles (Ny x Ne)
%           y: observation vector  (Ny x 1)
%          HC: matrix determining localization between obs- and model-space  (Ny x Nx)
%         HCH: matrix determining localization between obs- and obs-space  (Ny x Nx)
%        Neff: effective ensemble size for adaptive regularization coefficients
%       alpha: mixing parameter
%       var_y: obs error variance
%   kddm_flag: flag for performing probability mapping step (using Gaussian kernels) 
%
%  OUTPUT 
%        xmpf: posterior mean      (Nx x 1)
%           x: posterior particles (Nx x Ne)
%      e_flag: error flag

function [xmpf,x,e_flag] = pf_update(x,hx,y,HC,HCH,Neff,min_res,alpha,var_y,kddm_flag,qcpass,maxiter)

% Nothing to do if no obs pass QC
xmpf = mean(x')';
e_flag = 0;
if sum(qcpass) == length(y)
  return
end

% Get array dimensions
[Nx,Ne] = size(x);

% Remove obs that don't pass QC
ind = find(qcpass==1);
y(ind) = [];
hx(ind,:) = [];
HC(ind,:) = [];
HCH(ind,:) = [];
HCH(:,ind) = [];
Ny = length(y);

% Set parameter values for iterations
max_res = 1;         % current max residual

% Initialize variables for iterations
beta = ones(Nx,1);   % initial tempering coefficients
beta_y = ones(Ny,1); % initial tempering coefficients
beta_max = 1e100;     % maximum allowed tempering coefficient

res = ones(Nx,1);    % initial residuals
res_y = ones(Ny,1);  % initial residuals
niter = 0;           % initial iteration count
pf_infl = ones(Ny,1);
res_infl = ones(Ny,1);

% Remove minimum value
res = res - min_res;
res_y = res_y - min_res;

while max_res > 0 & min_res < 1

  niter = niter + 1;

  % Store copy of original particles
  xo = x;
  hxo = hx;

  % Initialize weighting matrices
  lomega = zeros(Nx,Ne);
  lomega_y = zeros(Ny,Ne);
  omega = ones(Nx,Ne)/Ne;
  omega_y = ones(Ny,Ne)/Ne;

  % Obs space weights
  for i = 1:Ny % Observation loop

    % Squared innovations
    d = (y(i) - hxo(i,:)).^2./(2.*var_y);
    d = d - min(d);

    % Obs space weights need to be capped to avoid collapse prior 
    % to calculating regularization coefficients 
    pf_infl = max(1,max(d)/200);
    pf_infl = find_obs_infl(d/pf_infl)*pf_infl;

    % Weight calculations
    d = d/pf_infl;
    wo(i,:) = exp( -d );
    wo(i,:) = wo(i,:) ./ sum(wo(i,:));

  end

  % Check for nans
  if sum(sum(isnan(wo))) > 0
    xmpf = ones(Nx,1)*nan;
    e_flag = 1; return  
  end

  % Calculate regularization coefficients
  [beta_y,res_y] = get_reg(Ny,Ny,Ne,HCH,wo,Neff,res_y,beta_max);
  [beta,res] = get_reg(Nx,Ny,Ne,HC,wo,Neff,res,beta_max);

  for i = 1:Ny % Observation loop

    % Skip if impact is low
    if (1 > 0.98*Ne*sum(wo(i,:).^2) ),continue, end
    wt = Ne*wo(i,:) - 1;

    % Model-space localized weighting vectors
    loc = HC(i,:);

    for j = 1:Nx

      if beta(j) == beta_max, continue, end

      dum = wt.*loc(j);
      ind = find(abs(dum)>0.3);
      dum(ind) = log( dum(ind) + 1 );

      lomega(j,:) = lomega(j,:) - dum;
      lomega(j,:) = lomega(j,:) - min(lomega(j,:));
    
    end

    % Obs-space localized weighting vectors
    loc = HCH(i,:);
    for j = 1:Ny

      if beta_y(j) == beta_max, continue, end

      dum = wt.*loc(j);
      ind = find(abs(dum)>0.3);
      dum(ind) = log( dum(ind) + 1 );

      lomega_y(j,:) = lomega_y(j,:) - dum;
      lomega_y(j,:) = lomega_y(j,:) - min(lomega_y(j,:));

    end

    % Normalize weights
    for n = 1:Ne
      omega(:,n) = exp(-lomega(:,n)./beta);
      omega_y(:,n) = exp(-lomega_y(:,n)./beta_y);
    end

    omegas_y = sum(omega_y')';
    omegas = sum(omega')';

    % Get posterior mean
    xmpf = 0;
    hxmpf = 0;

    for n = 1:Ne
      omega(:,n) = omega(:,n)./omegas;
      xmpf = xmpf + omega(:,n).*xo(:,n);
      omega_y(:,n) = omega_y(:,n)./omegas_y;
      hxmpf = hxmpf + omega_y(:,n).*hxo(:,n);
    end

    % Skip update step if few particles removed
    w = omega_y(i,:);
    if (1 > 0.98*Ne*sum(w.^2) ), continue, end

    % Calculate posterior variance using vector weights and original prior particles
    var_a = 0;
    var_a_y = 0;
    for n = 1:Ne
      var_a = var_a + omega(:,n).*(xo(:,n)-xmpf).^2;
      var_a_y = var_a_y + omega_y(:,n).*(hxo(:,n)-hxmpf).^2;
    end

    norm = ( 1 - sum(omega'.^2)' );
    var_a = var_a./norm;
    norm = ( 1 - sum(omega_y'.^2)' );
    var_a_y = var_a_y./norm;

    % Error check
    if sum(isnan(xmpf)) > 0
      e_flag = 1; return
    else
      e_flag = 0;
    end

    % Sample from original prior
    ind = sampling(hxo(i,:),omega_y(i,:),Ne);

    % Merge steps
    x  = pf_merge(x,xo(:,ind),HC(i,:)',Ne,xmpf,var_a,alpha);
    hx = pf_merge(hx,hxo(:,ind),HCH(i,:)',Ne,hxmpf,var_a_y,alpha);

  end

  % Use KDDM to update posterior samples
  if kddm_flag

    for j = 1:Nx
      x(j,:) = kddm(x(j,:),xo(j,:),omega(j,:));
    end

    xmpf = sum(x')'/Ne;

    for j = 1:Ny
      hx(j,:) = kddm(hx(j,:),hxo(j,:),omega_y(j,:));
    end

  end

  % Get max of residuals for all states
  max_res = max(res);
  if niter == maxiter, break, end

end % Iterations

% EnKF step (for hybrid only)
var_infl = ones(Ny,1)./min_res;

if max(min_res) > 0
  [xmpf,x] = enkf_update_tempered(x,hx,y,var_y,HC,HCH,0.95,var_infl);
end

