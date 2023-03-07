% Function that performs pf sampling step
%
% Input:
%        w <- weights used for determining sampling
%       Ne <- ensemble size

function ind = sampling(x,w,Ne);

%w = w.^0.5;  % sampling correction
%w = w./sum(w);

% Sort sample
[a,b] = sort(x);

% Apply deterministic sampling by taking value at every 1/Ne quantile
cum_weight = [0,cumsum(w(b))];

offset = 0.0;
base = 1/(Ne - offset)/2;

clear ind
k = 2;
for n = 1:Ne

  frac = base + (n - 1)/(Ne - offset);

  flag = 0;
  while flag == 0
    if ( cum_weight(k-1) < frac ) && ( frac <= cum_weight(k) )
      ind(n) = k-1;
      flag = 1;
    else
      k = k + 1;
    end
  end
end

% Unsort indices
ind = b(ind);

% Replace removed particles with duplicated particles
ind2 = ind*0;
for n = 1:Ne
  if sum(ind==n) ~= 0
    ind2(n) = n;
    dum = find(ind==n);
    ind(dum(1)) = [];
  end
end

ind0 = find(ind2==0);
ind2(ind0) = ind; ind = ind2;  

return

disp('------------')
for n = 1:Ne
%  if x(n) ~= x(ind(n))
    disp(['Replacing ',num2str(x(n)),' with ',num2str(x(ind(n)))])
%  end
end
