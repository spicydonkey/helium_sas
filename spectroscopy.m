%% Figure properties presets
figname = 'he_sas';
fontsize = 9;
linewidth = 1;
markersize = 5;

paperunits = 'centimeters';
papersize = [8 6];
paperposition = [0 0 papersize];


load('data\sas_data.mat');      

%% Plot scaled spectroscopy
fig = figure();

plot(f/1e12,sas_full,'linewidth',linewidth);
hold on;
%plot(f/1e12,sas_smooth+1,'linewidth',1.5);

axis tight;

x_lim = get(gca,'XLim');
y_lim = get(gca,'YLim');

set(gca,'XLim',x_lim+0.1*diff(x_lim)*[-1 1]);
set(gca,'YLim',y_lim+0.1*diff(y_lim)*[-1 1]);
%set(gca,'YLim',[-3.8,1.3]);

set(gca,'XTick',276.73:2e-3:276.74);
set(gca,'XMinorTick','on');

%title('Saturated absorption spectroscopy of Helium with ECDL');
xlabel('Frequency (THz)');
ylabel('Absorption signal (V)');

%%% Annotation
anno_t{1} = text(f_P2*1e-12,-1.5,'$2^3S_1-2^3P_2$',...
    'HorizontalAlignment','center','FontSize',10);
anno_t{2} = text((f_P2+df_P1_P2)*1e-12,0.2,'$2^3S_1-2^3P_1$',...
    'HorizontalAlignment','center','FontSize',10);
f_xover = f_P2+0.5*df_P1_P2;
anno_t{3} = text(f_xover*1e-12,0.75,'x/o',...
    'HorizontalAlignment','center','FontSize',10);

% %%% Error signal
% ind_range = 250:300;    % small range around the P2 peak
% plot(f(ind_range)/1e12,sas_full(ind_range),'k');


%% Postprocess
set(gca,'FontSize',fontsize);

fig.Units = paperunits;
fig.Position = paperposition;

fig.PaperSize = papersize;
fig.PaperUnits = paperunits;
fig.PaperPosition = paperposition;

%saveas(fig, [figname, '.eps'], 'psc2');     % save fig in cd


%% Zoom into 3P2 cooling transition
df_M = 1e-6 * (f - f_P2);       % diff freq from transition (MHz)
df_lim = 150;                   % half-width limit about 0 (MHz)

[~,I_min] = min(abs(df_M + df_lim));
[~,I_max] = min(abs(df_M - df_lim));

df_M_zoom = df_M(I_min:I_max);
y = sas_full(I_min:I_max);

% line profiles
f_gauss = @(b,x) b(4) + b(3) * exp(-(x-b(1)).^2 /(2*b(2)^2) );  % gaussian
f_lorentz = @(b,x) b(4) + b(3)./(1 + ((x - b(1))./b(2)).^2 );   % lorentzian

p0 = [0, 3, 1, -3];     % [x peak, width, amp, offset]

%%% FIT
fit_G = fitnlm(df_M_zoom,y,f_gauss,p0);
fit_L = fitnlm(df_M_zoom,y,f_lorentz,p0);

% fit prediction
xx = df_lim*linspace(-1,1,1e3);

% config
pred_nsig = 1;                  % n-sdev range: mu ± n*sigma
pred_conflvl = erf(pred_nsig/sqrt(2));  % confidence level
pred_alpha = 1-pred_conflvl;    % alpha: 100(1 – alpha)%

[yy_G, yy_G_ci] = predict(fit_G,xx','Alpha',pred_alpha,'Simultaneous',true);
[yy_L, yy_L_ci] = predict(fit_L,xx','Alpha',pred_alpha,'Simultaneous',true);

% fit summary
pfit_G = [fit_G.Coefficients.Estimate, fit_G.Coefficients.SE];
pfit_L = [fit_L.Coefficients.Estimate, fit_L.Coefficients.SE];

str_G = sprintf('$\\sigma =$ %0.2g $\\pm$ %0.1g MHz', pfit_G(2,1), pfit_G(2,2));
str_L = sprintf('$\\gamma =$ %0.2g $\\pm$ %0.1g MHz', pfit_L(2,1), pfit_L(2,2));

%%% plot
H = figure('Name','23P2_sas');
plot(df_M_zoom,y,'k*');

hold on;

p_fitG = plot(xx,yy_G,'b-','DisplayName',str_G);
p_fitG_ci = plot(xx,yy_G_ci,'b--');       % confidence interval

% p_fitL = plot(xx,yy_L,'r--','DisplayName',str_L);
% p_fitL_ci = plot(xx,yy_L_ci,'r--');       % confidence interval


xlabel('$\Delta f$ (MHz)');
ylabel('SAS signal (arb unit)');

xlim([-150,150]);

% legend
% legend([p_fitG,p_fitL]);
legend([p_fitG]);