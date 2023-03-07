% Function for merging prior and sampled particles during update steps
%
%  INPUT 

function xa = pf_merge(x,xs,loc,Ne,xmpf,var_a,alpha)

% coefficient for merging step
c = (1-loc)./loc;

% Subtract weight-estimated posterior mean from particles
xs = xs-xmpf;
x = x-xmpf;

% Calculate r1 and r2 coefficients for weight update equation
v1 = sum( xs'.^2 )';
v2 = sum( x'.^2 )';
v3 = sum( x'.*xs' )';

c2 = c.*c;
r1 = v1 + c2.*v2 + 2.*c.*v3;
r2 = c2./r1;

r1 = alpha*sqrt((Ne-1).*var_a./r1);
r2 = sqrt((Ne-1).*var_a.*r2);

% Calculate alpha2 which is added to r2 to fit posterior variance
m1 = sum( xs' )'/Ne;
m2 = sum( x' )'/Ne;
v1 = v1 - Ne.*m1.^2;
v2 = v2 - Ne.*m2.^2;
v3 = v3 - Ne.*m1.*m2;

T1 = v2;
T2 = 2.*( r1.*v3 + r2.*v2 );
T3 = v1.*r1.^2 + v2.*r2.^2 + 2.*v3.*r1.*r2 - (Ne-1).*var_a;
alpha2 = ( - T2 + sqrt( T2.^2 - 4.*T1.*T3 ) ) ./ (2.*T1);

r2 = r2 + alpha2;


% Generate localized posterior particles and calculate sample mean
pfm = 0;

for n = 1:Ne
  xa(:,n) = xmpf + r1.*xs(:,n) + r2.*x(:,n);
  pfm = pfm + xa(:,n);
end
pfm = pfm/Ne;

% Re-center particles on posterior mean
for n = 1:Ne
  xa(:,n) = xmpf + (xa(:,n) - pfm);
  % Skip update if nans exist (caused by zero variance)
  nanind = find(isnan(xa(:,n)));
  xa(nanind,n) = xmpf(nanind) + xs(nanind,n);
end
