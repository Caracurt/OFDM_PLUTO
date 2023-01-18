if exist('plutoTx')
    release(plutoTx)
end

if exist('plutoRx')
    release(plutoRx)
end

close all; clear all; clc;

%global plutoRadio
%plutoRadio = findPlutoRadio;

center_frequency = 1200e6; % 915e6 1.2e9
%sample_rate = 940e3; % 300e3
sample_rate = 940 * 1e3; % 300e3

N_length = 1024;
freq = 3/N_length;
freq1 = 256/N_length;
x_arr = 0:N_length-1;
sig_tx = 1*exp(1j*2*pi*freq*x_arr.') + 1*exp(1j*2*pi*freq1*x_arr.');
guard = zeros(N_length, 1);

IQmessage = [sig_tx; guard; sig_tx];
IQmessage = complex(IQmessage);

tx_signal = IQmessage;
plutoTx = plutoTxInit(center_frequency,sample_rate);
transmitRepeat(plutoTx,tx_signal);
status = 1;
pause(0.5)

samples_per_frame = 10 * ( 3*N_length );

%RTLSDRrx = RTLSDRRxInit(center_frequency,sample_rate,samples_per_frame);
%rx_signal = RTLSDRrx();

plutoRx = plutoRxInit(center_frequency,sample_rate,samples_per_frame);
rx_signal = plutoRx();


figure(1);
plot(10*log10(abs(rx_signal)))
figure(2);
plot((real(rx_signal)))
figure(3);
plot(abs(fft(rx_signal)))

release(plutoRx);



