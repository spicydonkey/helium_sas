%% Figure properties presets
figname = 'he_sas_mhf';
fontsize = 9;
fontsize_2 = 11;
linewidth = 1;
markersize = 20;

paperunits = 'centimeters';
% papersize = [8 6];
papersize = [17.6 6];
paperposition = [0 0 papersize];


load('data\sas_data.mat');      


%% fit stuff
pred_nsig = 1;                  % n-sdev range: mu ± n*sigma
pred_conflvl = erf(pred_nsig/sqrt(2));  % confidence level
pred_alpha = 1-pred_conflvl;    % alpha: 100(1 – alpha)%


%% Plot scaled spectroscopy
H_sas_mhf = figure('Name',figname);

p_data = plot(f/1e12,sas_full,'b-','linewidth',linewidth,'DisplayName','Data');
% plot(f/1e12,sas_full,'.','linewidth',linewidth);
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
ylabel('SAS signal (arb unit)');



% %%% Error signal
% ind_range = 250:300;    % small range around the P2 peak
% plot(f(ind_range)/1e12,sas_full(ind_range),'k');


%% fit Doppler broadened spectrum
f2x = @(f) (f - f_P2)*1e-9;
x2f = @(x) x*1e9 + f_P2;

% x = (f - f_P2)*1e-9;
x = f2x(f);
y = sas_full;


% absorption
model_p1p2 = @(b,x)   b(3)*exp(-(x-b(1)).^2/(2*b(2)^2)) ...
                    + b(6)*exp(-(x-b(4)).^2/(2*b(5)^2)) + b(7);     % sum of 2 gaussians (Doppler)
p0_p1p2 = [0, 1, -5, df_P1_P2*1e-9, 1, -4, 3];

fit_p1p2 = fitnlm(x,y,model_p1p2,p0_p1p2);


% abs 2 -- fix resonance freqs
model_p1p2_2 = @(b,x)   b(2)*exp(-(x-0).^2/(2*b(1)^2)) ...
                    + b(4)*exp(-(x-df_P1_P2*1e-9).^2/(2*b(3)^2)) + b(5);     % sum of 2 gaussians (Doppler)
p0_p1p2_2 = [1, -5, 1, -3, 3];

fit_p1p2_2 = fitnlm(x,y,model_p1p2_2,p0_p1p2_2);

% absorption 3 -- linear background
model_p1p2_3 = @(b,x)   b(3)*exp(-(x-b(1)).^2/(2*b(2)^2)) ...
                    + b(6)*exp(-(x-b(4)).^2/(2*b(5)^2)) + b(7) ...  % sum of 2 gaussians (Doppler)
                    + b(8)*x;       % linear background
p0_p1p2_3 = [0, 1, -5, df_P1_P2*1e-9, 1, -4, 3, 0.1];

fit_p1p2_3 = fitnlm(x,y,model_p1p2_3,p0_p1p2_3);


% saturated absorption
model_sas = @(b,x)  b(3)*exp(-(x-b(1)).^2/(2*b(2)^2)) ...
                  + b(6)*exp(-(x-b(4)).^2/(2*b(5)^2)) ...
                  + b(7)*exp(-(x-0.003774).^2/(2*0.006678^2)) ...  % SAS/xover peaks
                  + b(8)*exp(-(x-1.1615).^2/(2*0.003^2)) ...
                  + b(9)*exp(-(x-df_P1_P2*1e-9).^2/(2*0.003^2)) + b(10);     % sum of 2 gaussians
p0_sas = [0, 1, -5, df_P1_P2*1e-9, 1, -3, ...
    0.5, 0.2, 0.1, 5.8662];

fit_sas = fitnlm(x,y,model_sas,p0_sas);



% fit predictions
xfit = linspace_lim(min_max(x),1e3)';

[yfit, yfit_ci] = predict(fit_p1p2,xfit,'Alpha',pred_alpha,'Simultaneous',true);
[yfit2, yfit_ci2] = predict(fit_p1p2_2,xfit,'Alpha',pred_alpha,'Simultaneous',true);
[yfit3, ~] = predict(fit_p1p2_3,xfit,'Alpha',pred_alpha,'Simultaneous',true);
[yfit_sas, yfit_sas_ci] = predict(fit_sas,xfit,'Alpha',pred_alpha,'Simultaneous',true);


% plot
h_fits = figure('Name','fits');
hold on;
plot(x,y,'.','DisplayName','data');
plot(xfit,yfit,'-','DisplayName','Doppler');
plot(xfit,yfit2,'--','DisplayName','Doppler2');
plot(xfit,yfit3,':','DisplayName','Doppler3');
plot(xfit,yfit_sas,'-','DisplayName','SAS');
legend();


%% Main figure
figure(H_sas_mhf);

% plot a fitted one
pfit = plot(x2f(xfit)/1e12,yfit,'g-','linewidth',2*linewidth,...
    'DisplayName','Doppler fit');
uistack(pfit,'bottom');

set(gca,'FontSize',fontsize);

H_sas_mhf.Units = paperunits;
H_sas_mhf.Position = paperposition;

H_sas_mhf.PaperSize = papersize;
H_sas_mhf.PaperUnits = paperunits;
H_sas_mhf.PaperPosition = paperposition;



%%% Annotation
% arrows
f_anno = [f_P2,(f_P2+df_P1_P2),f_P2+0.5*df_P1_P2];
y_f = feval(fit_p1p2,f2x(f_anno));
y_i = y_f + 1.5;

text_anno = {'$2\,^3\textrm{S}_1 \rightarrow 2\,^3\textrm{P}_2$','$2\,^3\textrm{S}_1 \rightarrow 2\,^3\textrm{P}_1$','crossover'};

for ii=1:length(f_anno)
    text(f_anno(ii)/1e12,y_i(ii),text_anno{ii},...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',fontsize);
    arrow3([f_anno(ii)/1e12,y_i(ii)],[f_anno(ii)/1e12,y_f(ii)+0.5],'',1/2);
end
axis normal;


legend([p_data,pfit],'Location','southeast');

%saveas(H_sas_mhf, [figname, '.eps'], 'psc2');     % save H_sas_mhf in cd


%% Zoom into 3P2 cooling transition
df_M = 1e-6 * (f - f_P2);       % diff freq from transition (MHz)
df_lim = 150;                   % half-width limit about 0 (MHz)


% recenter zero by fit
df_0 = 3.9629;      % obtained from fit
df_M = df_M - df_0; 


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

[yy_G, yy_G_ci] = predict(fit_G,xx','Alpha',pred_alpha,'Simultaneous',true);
[yy_L, yy_L_ci] = predict(fit_L,xx','Alpha',pred_alpha,'Simultaneous',true);

% fit summary
pfit_G = [fit_G.Coefficients.Estimate, fit_G.Coefficients.SE];
pfit_L = [fit_L.Coefficients.Estimate, fit_L.Coefficients.SE];

str_G = sprintf('$\\sigma =$ %0.2g $\\pm$ %0.1g MHz', pfit_G(2,1), pfit_G(2,2));
str_L = sprintf('$\\gamma =$ %0.2g $\\pm$ %0.1g MHz', pfit_L(2,1), pfit_L(2,2));

str_G_lwidth = sprintf('%0.2g $\\pm$ %0.1g MHz', 2.3548*pfit_G(2,1), 2.3548*pfit_G(2,2));

%%% plot
H = figure('Name','23P2_sas');
plot(df_M_zoom,y,'b.','MarkerSize',markersize);

hold on;

p_fitG = plot(xx,yy_G,'k-','DisplayName','Gaussian fit');
p_fitG_ci = plot(xx,yy_G_ci,'k--');       % confidence interval

% p_fitL = plot(xx,yy_L,'r--','DisplayName',str_L);
% p_fitL_ci = plot(xx,yy_L_ci,'r--');       % confidence interval


xlabel('$\Delta f$ (MHz)');
ylabel('SAS signal (arb unit)');

xlim([-150,150]);


%%% Annotations
% FWHM linewidth
x_hmax = pfit_G(1,1) + sqrt(2*log(2)) * pfit_G(2,1) * [-1,1];     % x at half-maximum
y_hmax = pfit_G(4,1) + 0.5*pfit_G(3,1);   % y half maximum
larrow = 50;        % arrow length


pbaspect('manual');     % otherwise aspect ratio rescales
arrow3([x_hmax(1)-larrow,y_hmax],[x_hmax(1),y_hmax],'',2/3);
arrow3([x_hmax(2)+larrow,y_hmax],[x_hmax(2),y_hmax],'',2/3);
text(x_hmax(2)+larrow,y_hmax,str_G_lwidth,'FontSize',fontsize_2);


%%% legend
% legend([p_fitG,p_fitL]);
lgd = legend([p_fitG],'FontSize',fontsize_2);
set(gca,'FontSize',fontsize_2);
