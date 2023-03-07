% ==========================================================
% Main script for comparing the local particle filter of 
% Poterjoy (2022) and an ensemble square root filter 
% (Whitaker and Hamill 2002)
%
% This version supports three models: Lorenz (1963),
% Lorenz (1996), and two-scale Lorenz (1996)
%
% ==========================================================

% Link matlab libraries needed for experiments
da_libs

close all; clear all; warning off

%  ........................................................
% Choose a model with m_flag:   1 <== Lorenz 63
%				2 <== Lorenz 96
%				3 <== Lorenz 05 (model III)
m_flag = 2;

%  ........................................................
%  Specify filter with f_flag{}
%  ........................................................
%  :      f_flag options                                  :
%  ........................................................
%  : 1 <== EnKF	 (Whitaker and Hamill 2022)
%  : 2 <== Local PF (Poterjoy 2022)

f_flag{1} = 1;
f_flag{2} = 2;

% Experiment names
e_name{1} = 'EnKF'; e_name{2} = 'local PF';

%  --------------------------------------------------------
% | The first part of the code contains all the experiment |
% | parameters, such as model and observation information  |
% | and filter parameters.                                 |
%  --------------------------------------------------------

% --- Model parameters ---
T   = 1000;   % number of obs times
Ne  = 80;    % number of particles

% Experiment flags
plot_flag = 0;  % Set to 1 to plot prior values each filter time
disp_flag = 1;  % Set to 1 to display mean RMSEs each filter time
h_flag    = 0;  % Set to 0 for H(x) = x
                % Set to 1 for H(x) = x^2
                % Set to 2 for H(x) = log(|x|)

% --- Observation parameters ---
sig_y = 1;    % observation error
tau   = 1;    % model time steps between observation times
obf   = 4;    % observation spatial frequency: spacing between variables
obb   = 0;    % observation buffer: number of variables to skip when generating obs
              % - setting a non-zero obb will create a data void in the domain

% --- EnKF parameters ---
roi_kf   = 10;   % EnKF localization radius of influence 
inf_flag = 2;    % 0 ==> no inflation
                 % 1 ==> State-space adaptive inflation 
                 % 2 ==> RTPS with parameter gamma 
gamma    = 0.2;  % RTPS parameter (only relavent when inf_flag == 2)

% --- PF parameters ---
roi_pf = 10;         % PF localization radius of influence
alpha = 0.50;        % PF mixing coefficient
Neff =  0.40*Ne;     % effective ensemble size used for regularization
min_res = 0.0;       % minimum residual left after tempering with regularization
pf_kddm = 0;         % PF KDDM option

% Number of filter options
Nf = length(f_flag);

% Model parameters are determined by case block for selected model 
switch m_flag

  % Lorenz 63 model parameters
  case 1
    Nx  = 3;    % number of variables
    dt  = 0.01; % time step
    s = 10; 
    r = 28; 
    b = 8/3;

  % Lorenz 96 model parameters
  case 2, 
    Nx  = 40;   % number of variables
    dt  = 0.05; % time step
    F   = 8;    % forcing term (truth)
    Fe  = 8;    % forcing term (model)

  % Two-scale Lorenz 05 model parameters
  case 3, 

    Nx = 480;    % model variables
    dt = 0.05; % time step
    F  = 15;     % forcing term (truth)
    Fe = 15;    % forcing term (model)
    K  = 32;    % spatial smoothing parameter
    Im = 12;
    b  = 10.0;
    c  = 2.5;

    % Adjust localization to be consistent with typical L96 model
    roi_kf = roi_kf*24;
    roi_pf = roi_pf*24;

end

% Use same random numbers each experiment
rng(1); 

%  --------------------------------------------------------
% | The next part of the code sets up the experiment. It   |
% | generates the measurement operator, localization       |
% | matrix, and obs error covariance matrix, spins up the  |
% | model, creates a truth simulation and observations,    |
% | and generates the initial prior ensemble. The truth    |
% | simulation is formed by running the model 1000 time    |
% | steps + the length of the experiment specified with T. |
% | The first prior ensemble comes from an ensemble        | 
% | forecast initialized from the randomly perturbed truth |
% | state.                                                 |
%  --------------------------------------------------------

% Generate observation error covariance matrix
var_y  = sig_y^2; 

numobs = ceil((Nx-2*obb)/obf);
R      = eye(numobs)*var_y;
R_i    = inv(R);

% Define H
H  = eye(Nx); 
H  = H(obb+1:obf:Nx-obb,:);
Ny = length(H(:,1));


% Correlation matrix for localization
C_kf  = gen_be_periodic(roi_kf,1,Nx,1);
C_pf  = gen_be_periodic(roi_pf,1,Nx,1);

% Apply interpolation part of measurement operator 
C_kf = H*C_kf;
C_pf = H*C_pf;

% Define domain
xd = [1:Nx]';

% Initialize model for spinup period
xt(1:Nx,1) = 3*sin([1:Nx]/(6*2*pi));

% Spin up initial truth state
switch m_flag
  case 1
    xt = M_nl_l63(xt,dt,1000,s,r,b);
  case 2
    xt = M_nl_l96(xt,dt,1000,F);
  case 3
    xt = M_nl_l05III(xt,dt,1000*0.05/dt,K,Im,b,c,F);
end

% Start parallel run
%delete(gcp('nocreate'))
%poolobj = parpool(10);

% Run initial ensemble forecast
%parfor n = 1:Ne
for n = 1:Ne

  dum = xt + 1*randn(Nx,1);
  switch m_flag
    case 1
      xi(:,n) = M_nl_l63(dum,dt,100,s,r,b);
    case 2
      xi(:,n) = M_nl_l96(dum,dt,100,Fe);
    case 3
      xi(:,n) = M_nl_l05III(dum,dt,50,K,Im,b,c,F);
  end

end

% Generate Truth
switch m_flag

  case 1
    xt = M_nl_l63(xt,dt,100,s,r,b);
    for t = 2:T
      xt(:,t) = M_nl_l63(xt(:,t-1),dt,tau,s,r,b);
    end

  case 2
    xt = M_nl_l96(xt,dt,100,F);
    for t = 2:T
      xt(:,t) = M_nl_l96(xt(:,t-1),dt,tau,F);
    end

  case 3
    xt = M_nl_l05III(xt,dt,100,K,Im,b,c,F);
    for t = 2:T
      xt(:,t) = M_nl_l05III(xt(:,t-1),dt,tau,K,Im,b,c,F);
    end

end

% Create synthetic obs from truth and random errors
dum = randn(T,Nx)'*sig_y;

switch h_flag
  case 0 % ---   H(x) = x + eps   ---
    Y = H*( xt + dum );
  case 1 % ---   H(x) = x^2 + eps   ---
    Y = H*( (xt.^2 + dum ) );  
  case 2 % ---   H(x) = log(abs(x)) + eps   ---
    Y = H*( log(abs(xt + dum )) );
end

% Initialize prior ensemble for all experiments
for f = 1:Nf
  x{f} = xi;
end

% Initialize inflation values for adaptive inflation
for f = 1:Nf
  prior_inf{f} = ones(Nx,1);
  prior_inf_y{f} = ones(Ny,1);
  var_inf = 0.8;
end

%  --------------------------------------------------------
% | This part of the code loops through the observation    |
% | times. It applies a filter to update the ensemble      |
% | then runs an ensemble forecast from updated members.   |
%  --------------------------------------------------------

e_flag = 0;
for t = 1:T % Time loop

  % Plot prior information for each filter time
  if plot_flag

    for f = 1:Nf
      figure(f)
      % Plot prior particles 
      subplot(3,1,1); hold off
      for n = 1:Ne
        plot(x{f}(:,n),'color','b','linewidth',2); hold on;
      end      

      title(['Prior particles (',e_name{f},')'],'fontsize',20)

      % Plot observations
      scatter(H*xd,Y(:,t),'k','linewidth',2);

      % Plot truth
      plot(xt(:,t),'color','g','linewidth',2);
      xlim([0,Nx+1]);

      % Plot histogram for variable 1
      subplot(3,1,2); hold off;
      hist(x{f}(1,:),Nx)
      title('X1 prior histogram','fontsize',20)
    end

  end

  % Loop through experiments and perform update based on specified options
  for f = 1:Nf

    % --------------------------
    % ------ Filter Step -------
    % --------------------------

    % Obs-space priors
    switch h_flag
      case 0; hx = H*x{f};
      case 1; hx = H*( (x{f}.^2 ) );
      case 2; hx = H*log(abs(x{f}));
    end

    % QC step
    qcpass = zeros(1,Ny);
    for i = 1:Ny
      d = abs(Y(i,t) - mean(hx(i,:)));
      if d > 4 * sqrt( var(hx(i,:)) + var_y  )
        qcpass(i) = 1;
      end
    end
    clear d

    % Call filter 
    switch f_flag{f}

      case 1 % EnKF update step

        [xm{f},x{f},e_flag,prior_inf{f},prior_inf_y{f}] = enkf_update(x{f},hx,Y(:,t),var_y,C_kf,C_kf*H',inf_flag,prior_inf{f},prior_inf_y{f},var_inf,gamma,qcpass);

      case 2 % Local PF update step

        [xm{f},x{f},e_flag] = pf_update(x{f},hx,Y(:,t),C_pf,C_pf*H',Neff,min_res,alpha,var_y,pf_kddm,qcpass,3);

    end

    % Set all values to nan if filter fails
    if max(e_flag) == 1
      xm{f} = xm{f}*nan;
      x{f} = x{f}*nan;
    end

  end  % DA experiment loop

  % Plot posterior information for each filter time
  if plot_flag

    for f = 1:Nf

      figure(f)
      subplot(3,1,3); hold off
      for n = 1:Ne
        plot(x{f}(:,n),'color','b','linewidth',2); hold on;
        title('Posterior particles','fontsize',20)
        xlim([0,Nx+1]);
      end
      scatter(H*xd,Y(:,t),'k','linewidth',2);
      plot(xt(:,t),'color','g','linewidth',2);
      plot(xm{f},'color','r','linestyle','--','linewidth',2);
      pause(1)

    end
  end

%  --------------------------------------------------------
% | The next part of code in the time loop calculates and  |
% | saves some basic statistics from the experiments.      |
%  --------------------------------------------------------

  % Keep track of posterior errors during experiments

  if disp_flag
    tstep = sprintf('%2.2f',t);
    textl = ['Time = ',tstep];
  end

  for f = 1:Nf

    % RMSE errors
    dif = xt(:,t) - xm{f}; 
    err{f}(t) = sqrt(mean(dif.*dif));

    % Spread
    xvar = 0;
    for n = 1:Ne
      xvar = xvar + ((x{f}(:,n) - xm{f})).*((x{f}(:,n) - xm{f}))./(Ne-1);
    end

    sig{f}(t) = sqrt(mean(xvar));

    % Values for display
    if disp_flag
      rmse = sprintf('%2.4f',err{f}(t));
      textl = [textl,', RMSE (',e_name{f},') = ',rmse];
    end

    % Print errors and spread at each time
    if disp_flag && f == Nf
      disp(textl)
    end

  end

  % Run ensemble forecast for next cycle
  for f = 1:Nf
    dum = x{f};
%    parfor n = 1:Ne
    for n = 1:Ne
      switch m_flag
        case 1 % Lorenz 63
          dum(:,n) = M_nl_l63(dum(:,n),dt,tau,s,r,b);
        case 2 % Lorenz 96
          dum(:,n) = M_nl_l96(dum(:,n),dt,tau,Fe);
        case 3 % Two-scale Lorenz 05
          dum(:,n) = M_nl_l05III(dum(:,n),dt,tau,K,Im,b,c,F);
      end
    end
    x{f} = dum;
  end

end % Time loop

%  --------------------------------------------------------
% | The last part makes figures and saves data from each   |
% | experiment.                                            |
%  --------------------------------------------------------

fig = figure(2); hold off;

cl = colormap(lines(Nf));
for f = 1:Nf
  plot(err{f},'color',cl(f,:),'linewidth',2,'linestyle','-'); hold on;
  plot(sig{f},'color',cl(f,:),'linewidth',2,'linestyle','--');

  % Take average of rmse and spread
  mea(f)  = mean(err{f}(50:end));
  sigd(f) = mean(sig{f}(50:end));

  % Legend
  leglab{2*f-1} = [e_name{f},' RMSE (average: ',num2str(mea(f)),')'];
  leglab{2*f} = [e_name{f},' spread (average: ',num2str(sigd(f)),')'];

end

xlim([0,t+1]);
ylabel('Mean RMSE/spread','fontsize',20)
xlabel('Cycle number','fontsize',20)
set(gca,'fontsize',16);

lh = legend(leglab,'location','northwest');
set(lh,'box','off')

disp(' ')
for f = 1:Nf
  disp(['Time-average ',e_name{f},' RMSE / spread: ',num2str(mea(f)), ' / ',num2str(sigd(f))])
end
disp(' ')

return

% Save data
data.err = err;
data.sig = sig;
ofile = ['DATA/pf_rmse_sprd_sig_y_',num2str(sig_y),'.mat'];
fid = fopen(ofile,'w');
disp(['Saving data to ',ofile])
save(ofile,'data');
fclose(fid);

% Save plot
set(gca,'ylim',[0,5])
set(fig,'PaperSize',[12,5])
set(fig,'PaperPositionMode','manual'); 
set(fig,'PaperPosition',[0,0.1,12,5]);
ofile = ['FIGS/pf_rmse_sprd_',num2str(sig_y),'.pdf'];
disp(['Saving figure to ',ofile])
saveas(fig, ofile,'pdf')
