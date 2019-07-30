
clear all;


%% Figure properties presets
figname = 'he_sas';
fontsize = 9;
linewidth = 1;
markersize = 5;

paperunits = 'centimeters';
papersize = [8 6];
paperposition = [0 0 papersize];


load('C:\Users\David\Dropbox\PhD\lab\ecdl\He_spectroscopy\data\sas_data.mat');

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