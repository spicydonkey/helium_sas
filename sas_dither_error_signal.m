% LIA error signal from dither lock to SAS
% DKS
% 2019-07-30


%% Figure properties presets
% figname = 'he_sas_mhf';
fontsize = 11;
linewidth = 1;
markersize = 8;

paperunits = 'centimeters';
% papersize = [8 6];
% papersize = [17.6 6];
% paperposition = [0 0 papersize];




%% load data
data = read_osci_rigol('data/raw/exp_data/he_spec_5.csv');

t = data(:,1);
v = data(:,2);
sas = data(:,3);
err = data(:,4);


%% calibrate frequency
% post processed result
f0_fit = 8.5853e6;   % P2 location from fit (MHz)


% maually I find the 355th sample is the reference 2^3S_1 --> 2^3P_2 transition
I_P2 = 355;
f_P2 = 276.731672960526e12;       % transition in Hz

I_P1 = 800;
f_P1 = f_P2 + 2.29e9;

dfdv = (f_P1 - f_P2)/(v(I_P1) - v(I_P2));
f = (v - v(I_P2))*dfdv;      


dvdt = (v(I_P1) - v(I_P2))/(t(I_P1) - t(I_P2));
dfdt = dfdv * dvdt;
f = (t - t(I_P2))*dfdt;

f = f - f0_fit;     % shift such that zero-crossing is at x = 0

% filter single pass
I_filter = 209:957;

ff = f(I_filter);
sasf = sas(I_filter);
errf = err(I_filter);


%% plot
dy_arb = -0.3;      % arbitrary shift to y-axis to zero err signal at lock
errf0 = errf+dy_arb;


H = figure('Name','err_signal');

hold on;

% plot(ff/1e9,sasf,'b.');
plot(ff/1e9,errf0,'r.-');

xlabel('$\Delta f$ (GHz)');
ylabel('error signal (arb unit)');

box on;
xlim([min(ff),max(ff)]/1e9);
set(gca,'FontSize',fontsize);


%% Zoom-in plot
flim_inset = 0.2e9*[-1,1];
idxlim_inset = idxNearest(ff,flim_inset);       % get index of range

% get data for inset
f_inset = ff(idxlim_inset(1):idxlim_inset(2));
errf0_inset = errf0(idxlim_inset(1):idxlim_inset(2));


% Model: first derivative of gaussian
f_dgauss = @(b,x) - b(3) * (x - b(1)) .* exp( -(x - b(1)).^2 / (2*b(2)^2) );
p0 = [10,10,-1];        % initial param: in MHz units

% fit model to data
fit_dgauss = fitnlm(f_inset*1e-6,errf0_inset,f_dgauss,p0);


% fit prediction
% config
pred_nsig = 1;                  % n-sdev range: mu ± n*sigma
pred_conflvl = erf(pred_nsig/sqrt(2));  % confidence level
pred_alpha = 1-pred_conflvl;    % alpha: 100(1 – alpha)%

xx = 1e-6*linspace(flim_inset(1),flim_inset(2),1e3);        % in MHz
[yy, yy_ci] = predict(fit_dgauss,xx','Alpha',pred_alpha,'Simultaneous',true);


%%% plot
H=figure('Name','err_zoom_inset','Position',[0,0,325,154]);

hold on;

% FIT (GHz units)
plot(1e-3*xx,yy,'k-');       
% plot(xx,yy_ci,'r--');

% data
plot(1e-9*f_inset,errf0_inset,'r.','MarkerSize',markersize);        

box on;
axis tight;
set(gca,'FontSize',fontsize);