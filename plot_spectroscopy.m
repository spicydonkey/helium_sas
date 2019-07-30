filename = 'C:\Users\David\Dropbox\PhD\lab\ecdl\He_spectroscopy\data\raw\exp_data2\spec_opt_2.csv';

data = read_osci_rigol(filename);

figure();
plot(data(:,1),data(:,2));
hold on;
plot(data(:,1),data(:,3));
%plot(data(:,1),data(:,4));
grid on;

title('helium spectroscopy with ECDL');
xlabel('frequency (arb. units)');
ylabel('signal (arb. units)');

%legend('V_{PZT}','Doppler-free spectroscopy','LIA output');
legend('$V_{\textrm{PZT}}$','Doppler-free spectroscopy');