% Apply Kernal Density Distribution Mapping (adapted from code provided by Seth McGinnis)
% This version calculates quantiles of prior explicitely and inverts posterior cdf numerically
%
% Input:
%        x <- first guess sample
%       xo <- original prior sample
%        w <- weights calculated using original prior   

function xa = kddm_test(x,xo,w)

warning off;

Ne = length(w);

% Set bandwidth so the resulting kernel-estimated pdf has
% the same variance as the weighted estimate
%sig = sqrt(var(xo))*0.25;

sig = (max(x) - min(x))/6;
%sig = (max(x) - min(x))/12;

npoints = 300;

% Domain for defining posterior cdf
xmin = min(min(xo),min(x));
xmax = max(max(xo),max(x));

incr = (xmax-xmin)/(npoints-1);
xd = [xmin:incr:xmax];

qf = 0;
cdfxf = 0;
cdfxa = 0;

for n = 1:Ne

  % Approximate prior cdf
  %  cdfxf = cdfxf + (1/Ne)*( 1 + erf( (xd - x(n))/sqrt(2)/sig ) )/2;

  % Get quantiles of prior using sum of Gaussian cdfs evaluated at prior points
  qf = qf + ( 1 + erf( ( x - x(n) )/(sqrt(2)*sig) ) )/(2*Ne);
  cdfxa = cdfxa + w(n) * ( 1 + erf( (xd - xo(n))/(sqrt(2)*sig) ) )/2;

end

% Fix quantiles that fall outside domain xd
%if min(qf) < min(cdfxa)
%  r = [1:-1/(Ne-1):0];
%  [a,b] = sort(x);
%  qf(b) = qf(b) + r.*( min(cdfxa) - min(qf) );
%end
%if max(qf) > max(cdfxa)
%  r = [0:1/(Ne-1):1];
%  [a,b] = sort(x);
%  qf(b) = qf(b) + r.*( max(cdfxa) - max(qf) );
%end

% Invert cdfxa to get values at quantiles
xa = lininterp(cdfxa,xd,qf);


if var(xa) < 1e-8
  figure(3); hold off
%  plot(xd,cdfxf,'b'); hold on;
  plot(xd,cdfxa,'r'); hold on;
  scatter(x,qf,'b');
  scatter(xa,qf,'r');
  clear all; return
end

return

figure(4); hold off; plot(cdfxa); pause(3)


if sum(isnan(qf)) > 0
%if 1 == 1
  figure(3); hold off
%  plot(xd,cdfxf,'b'); hold on;
  plot(xd,cdfxa,'r');
  scatter(x,qf,'b');
  scatter(xa,qf,'r');

  for i = 1:Ne
    line([xo(i),xo(i)],[0,1])
  end

  figure(4); hold off;
  scatter(xo,ones(1,Ne)*0.2,'filled'); hold on
  [a,b] = sort(xo);
  scatter(xo(b),w(b));
  scatter(x,ones(1,Ne)*0.5,'filled');
  scatter(xa,ones(1,Ne)*0.8,'filled');

  set(gca,'ylim',[0,1.5])
%  clear all; return

end


% Scale to fit IS mean and variance
%xa = xm + ( xa - mean(xa) ).*sqrt(vm/var(xa));
