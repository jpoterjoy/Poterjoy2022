% Function for finding beta based on effective ensemble size

function beta = find_beta(sum_exp,Neff)

% Original weights
Ne = length(sum_exp);

beta_max = max(1,10*max(sum_exp));

w = exp(-sum_exp);
ws = sum(w);
if ws > 0
  w = w ./ ws;
  Neff_init = 1/sum(w.^2);
else
  Neff_init = 1;
end

% TEMP
%mix = 0.5;
%Neff = Neff_init + mix*(Ne - Neff_init);

if Neff == 1, beta = 1; return, end

if Neff_init < Neff || ws == 0

  % Initial start and end bounds
  ks = 1; ke = beta_max;

  % Apply bisection method to find k
  tol = 1e-5;

  for i = 1:1000

    % Evaluate function at end points
    w = exp( -sum_exp/ks );
    w = w./sum(w);    
    fks = Neff - 1 ./ sum(w.^2);
    if isnan(fks), fks = Neff-1; end

    w = exp( -sum_exp/ke );
    w = w./sum(w);
    fke = Neff - 1 ./ sum(w.^2);

    % Evaluate function at middle point
    km = (ke + ks) / 2;
    w = exp( -sum_exp/km );
    w = w./sum(w);
    fkm = Neff - 1 ./ sum(w.^2);
    if isnan(fkm), fkm = Neff-1; end

    %disp(['iteration: ',num2str(i),', km: ',num2str(km),', fks: ',num2str(fks),', fke: ',num2str(fke)])

    % Exit criteria
    if (ke-ks)/2 < tol, break, end

    if fkm*fks > 0
      ks = km;
    else
      ke = km;
    end

  end

  % Get beta from k
  beta = km;
  w = exp( -sum_exp/beta );
  w = w ./ sum(w);
  Nf = 1./sum(w.^2);

  %disp(['Target Neff: ',num2str(Neff),' | Actual Neff: ',num2str(Nf)])


  % In extreme cases, numerical errors can lead to the wrong result
  if ( (Nf <= Neff-1) || isnan(Nf) )
    disp(['WARNING! Neff is ',num2str(Nf),' but target is ',num2str(Neff)])
    beta = beta_max;
  % clear all; return
  end

%  disp(['number of iterations: ',num2str(i)])

else

  beta = 1;

end
