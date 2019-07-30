% LIA error signal from dither lock to SAS
% DKS
% 2019-07-30


%% load data
data = read_osci_rigol('data/raw/exp_data/he_spec_5.csv');

t = data(:,1);
v = data(:,2);
sas = data(:,3);
err = data(:,4);

%% calibrate frequency
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

% filter single pass
I_filter = 209:957;

ff = f(I_filter);
sasf = sas(I_filter);
errf = err(I_filter);


%% plot
dy_arb = -0.3;      % arbitrary shift to y-axis to zero err signal at lock
    
H = figure('Name','err_signal');

hold on;

% plot(ff/1e9,sasf,'b.');
plot(ff/1e9,errf+dy_arb,'r.-');

xlabel('$\Delta f$ (GHz)');
ylabel('error signal (arb unit)');

box on;
xlim([min(ff),max(ff)]/1e9);