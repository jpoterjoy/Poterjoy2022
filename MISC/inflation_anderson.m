% Anderson (2009) adaptive inflation: function modfied from Yue Ying (2014)
%
%  INPUT 
%           x: prior ensemble (Nx x Ne)
%          hx: obs-space prior ensemble (Ny x Ne)
%           y: observation vector  (Ny x 1)
%       var_y: obs error variance
%         inf: prior mean of inf coefficients (Nx x 1)
%     var_inf: prior variance of inf coefficients
%           C: matrix determining localization (Ny x Nx)
%        damp: damping parameter
%
%  OUTPUT 
%          xm: prior mean (Nx x 1)
%          xp: inflated prior ensemble perturbations (Nx x Ne)
%         inf: posterior mean of inf coefficients (Nx x 1)

function [xm,xp,inf]=inflation_anderson(x,hx,y,inf,var_inf,var_y,HC,damp,qcpass)

% Get array dimensions
[Nx,Ne] = size(x);
Ny = length(y);

% Separate mean and perturbations
xm = sum(x')'/Ne;
xp = x - xm;

for i = 1:Ny

  % QC check
  if qcpass(i) > 0, continue, end


  hxm = mean(hx(i,:));
  hxi = hx(i,:) - hxm;

  dist = abs(hxm - y(i));
  var = hxi*hxi'/(Ne-1);

  cov = xp*hxi'/(Ne-1);

  gama=HC(i,:).*cov/var;
  inf_o=(1+gama.*(sqrt(inf)-1)).^2;

  th=sqrt(inf_o*var+var_y);
  lm=exp(-0.5*(dist./th).^2)./(sqrt(2.*pi)*th);
  lp=lm.*((dist./th).^2-1)./th;
  lp=lp*0.5*var.*gama.*(1-gama+gama.*sqrt(inf))./(th.*sqrt(inf));
  bb=lm./lp-2*inf;
  cc=inf.^2-var_inf-lm.*inf./lp;
  inf1=0.5*(-bb+sqrt(bb.^2-4*cc));
  inf2=0.5*(-bb-sqrt(bb.^2-4*cc));
  for j=1:Nx
    if(gama(j)>0.0)
      if(abs(inf1(j)-inf(j)) > abs(inf2(j)-inf(j)))
        if(inf2(j)>0)
            inf(j)=inf2(j);
        end
      else
        if(inf1(j)>0)
          inf(j)=inf1(j);
        end
      end
    end
  end
end

% Limit inflation to values greater than 1
inf(find(inf<1)) = 1;

% Damping
inf = (inf - 1)*damp + 1;

% Inflate ensemble
for n = 1:Ne
  xp(:,n) = xp(:,n).*sqrt(inf);
end
