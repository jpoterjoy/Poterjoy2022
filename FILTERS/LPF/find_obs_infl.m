% Function for finding obs inflation based on minimum normalized weight

function alpha = find_obs_infl(lw)

% Original weights
Ne = length(lw);
minwt = 1e-10 / Ne;

% Find starting minimum weight
w = exp(-lw);
ws = sum(w);
if ws > 0
  w = w ./ ws;
  minw = min(w);
else
  minw = 0;
end

% Return alpha of 1 if threshold is already reached
if minw > minwt, alpha = 1; return, end

% Initial start and end bounds
ks = 1; ke = max(lw);

% Apply bisection method to find k
tol = 1e-4;

for i = 1:1000

  % Evaluate function at end points
  w = exp( -lw/ks );
  w = w./sum(w);    
  fks = min(w) - minwt;
  if isnan(fks), fks = -minwt; end

  w = exp( -lw/ke );
  w = w./sum(w);
  fke = min(w) - minwt;

  % Evaluate function at middle point
  km = (ke + ks) / 2;
  w = exp( -lw/km );
  w = w./sum(w);
  fkm = min(w) - minwt;

  %disp(['iteration: ',num2str(i),', km: ',num2str(km),', fks: ',num2str(fks),', fke: ',num2str(fke)])

  % Exit criteria
  if (ke-ks)/2 < tol, break, end

  if fkm*fks > 0
    ks = km;
  else
    ke = km;
  end

end

alpha = km;
w = exp( -lw/alpha );
w = w ./ sum(w);

%if isnan(min(w)), alpha = max(lw)/20; end

%disp(['Target min w: ',num2str(minwt),' | Actual min w: ',num2str(min(w))])
%disp(['number of iterations: ',num2str(i)])
