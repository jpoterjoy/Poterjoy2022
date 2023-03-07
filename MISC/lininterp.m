% Linear interpolation

function xi = lininterp(fx,x,fxi);

for n = 1:length(fxi)

  if fxi(n) >= fx(end)
    xi(n) = x(end);
    continue
  end

  if fxi(n) <= fx(1)
    xi(n) = x(1);
    continue
  end

  [dum,m] = min( fx < fxi(n) );

%    disp(['cda(m) | cda(m-1) | q ',num2str(fx(m)),' | ',num2str(fx(m-1)),' | ',num2str(fxi(n)),' | '])
%if fxi(n) > fx(m) || fxi(n) < fx(m-1), clear all; disp('ERROR'), return, end

  d = fx(m) - fx(m-1);

  if d < 1e-10
    xi(n) = x(m-1);
  else
    w1 = ( fx(m) - fxi(n) ) / d;
    w2 = ( fxi(n) -  fx(m-1) ) / d;
    xi(n) = w1*x(m-1) + w2*x(m);
  end

end
