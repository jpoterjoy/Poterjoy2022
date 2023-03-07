% Function for generating a univariate kernel estimate of a pdf
%
% Input:
%       xm <- sample
%        w <- sample weights used to define posterior pdf

function [fx,x] = kernel_density(xm,w)

w = w./sum(w);

mea = sum(w.*xm);
sig2 = sqrt( sum(w.*(xm-mea).^2) ./ ( 1 - sum(w.^2) ) );
%sig = 1;
%sig = sig2/2;
sig = sig2/10;
%sig = sig2/100;


% Parameters for online likelihood estimate (BORROWED FROM PF_UPDATE)
%mea = sum(w.*xm);
%mea2 = sum(w.*(xm.^2));
%sig2 = sum(w.*(xm-mea).^2) ./ ( 1 - sum(w.^2) );
%diff = sig2 - (mea2 - mea.^2);
%sig2 = sig2 - diff;
%sig = sig2/300;


% Domain for perfoming the mapping
%x = [ min(xm) - 4 : (8 + max(xm) - min(xm) )/1000 : max(xm) + 4 ];
%x = [ min(xm) - 4 : (8 + max(xm) - min(xm) )/5000 : max(xm) + 4 ];
x = xm;

xv = sqrt(var(xm));
%xmin = min(xm); xmax = max(xm);
%range = (xmax-xmin)*5/4;
%incr1 = 3*range/500;
%incr2 = 3*range/5000;
%x = [xmin-range:incr1:xmin,xmin+incr2:incr2:xmax,xmax+incr1:incr1:xmax+range];

%close all;
%scatter(x1,ones(1,length(x1)),'b'); hold on;
%scatter(x,0.5*ones(1,length(x)),'r');
%set(gca,'ylim',[0,2])
%return

% Bandwidth for Gaussian kernels is made somewhat adaptive by
% defining it as a function of distance between variables.
% This is done by setting the bandwith equal to the distance
% between the current variable and the Kth nearest variable.

% Number of points to consider for variable bw
K = floor(length(xm)/8);

fx = 0;
for i = 1:length(xm)

  % Distance between current data point and all other points
  dis = abs(xm(i) - xm);  

  % Sort distances
  dis = sort(dis);

  % Use average distance between Kth nearest point and current point for bw
%  sig = 0.1*mean(dis(1:K));
%sig = 4*sqrt( sum (dis(1:K) ) )/K;

  % Sum Gaussian kernels
  fx = fx + w(i)*exp(- 0.5 * (x-xm(i)).^2/sig^2 )/sqrt(2*pi)/sig;

end

