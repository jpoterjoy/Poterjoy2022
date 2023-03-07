% Function for merging prior and sampled particles during update steps
%
%  INPUT 

function [beta,res] = get_reg(Nx,Ny,Ne,C,hw,Neff,res,beta_max)

for j = 1:Nx

  % Nothing to do if min_res has already been reached
  if res(j) <= 0, beta(j,1) = beta_max; continue, end

  wo = 0;
  for i = 1:Ny
    dum = (Ne*hw(i,:) - 1).*C(i,j);
    ind = find(abs(dum)>0.3);
    dum(ind) = log( dum(ind) + 1 );
    wo = wo - dum;
    wo = wo - min(wo);
  end

  % Calculate beta from log of weights
  beta(j,1) = find_beta(wo,Neff);

  % Fix beta if value exceeds residual
  if res(j) < 1/beta(j)
    beta(j) = 1/res(j);
    res(j) = 0;
  else
    res(j) = res(j) - 1/beta(j);
  end

  beta(j) = min(beta(j),beta_max);

end

