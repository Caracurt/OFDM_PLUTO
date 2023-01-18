% transmit and receive harmonic function on sdr
if exist('plutoTx')
    release(plutoTx)
end

if exist('RTLSDRrx')
    release(RTLSDRrx)
end

close all; clear all; clc;

global plutoRadio

center_frequency = 1.2e9;
sample_rate = 940e3;
data_freq = sample_rate/10; % frequecny of sin
%data_rate = 

data_length = 128;
guard_length = 128;
freq_data = 0/data_length;
x_arr = 0:data_length-1;
x_arr = x_arr';
sig_tx = exp(1j*2*pi*x_arr*freq_data);
guard_int = zeros(guard_length, 1);
%tx_mess = complex( [sig_tx; guard_int; -sig_tx; guard_int] );
tx_mess = complex( [sig_tx; guard_int; -sig_tx; guard_int; guard_int; guard_int] );

%tx_mess = complex( [sig_tx; sig_tx; sig_tx; sig_tx] );

sig_len = length(tx_mess);

figure(10);
plot(real(tx_mess))

% Pluto transmit

% prepare message 


% upsampling and RRC
%sps = sample_rate/data_rate;

%txFilter = comm.RaisedCosineTransmitFilter( ...
%    'OutputSamplesPerSymbol',sps,'FilterSpanInSymbols',2);

%tx_signal = txFilter(IQmessage);

plutoRadio = findPlutoRadio;

plutoTx = plutoTxInit(center_frequency,sample_rate);

transmitRepeat(plutoTx,tx_mess);
% actual Tx

% pause
pause(1.0)

%return;

samples_per_frame = 20*sig_len;
RTLSDRrx = RTLSDRRxInit(center_frequency,sample_rate,samples_per_frame);

rx_signal = RTLSDRrx();
release(RTLSDRrx);
%release(plutoTx);

figure(20);
plot(10*log10(abs(rx_signal)));

figure(1);
title('Raw');
plot(real(rx_signal));

pos_high = find(abs(rx_signal) > 0.5);
pos_high = pos_high(1);
sig_cut = rx_signal(pos_high : pos_high + data_length/2 - 1);

data_len_fft = 8*length(sig_cut);
fft_sig = fft(sig_cut, data_len_fft);
figure(111);
plot(abs(fft_sig));

[~, pos_freq] = max(abs(fft_sig(1:end)));

freq_corr = -(pos_freq - 1)/data_len_fft;

time_arr = 0:length(rx_signal)-1;
time_arr = time_arr';
rx_signal_0 = rx_signal .* exp(1j*2*pi*freq_corr*time_arr);


figure(15);
title('Raw1');
plot(real(rx_signal_0));

sig_cut1 = rx_signal_0(pos_high : pos_high + data_length/2 - 1);

figure(16);
title('Raw1');
plot(imag(rx_signal_0));

if 0

% coarse sync
coarse = comm.CoarseFrequencyCompensator('SampleRate',sample_rate, ...
    'Modulation','BPSK','FrequencyResolution',1);
fine = comm.CarrierSynchronizer( ...
    'SamplesPerSymbol',1,'Modulation','BPSK');

[rx_coarse,estfreqoff] = coarse(rx_signal);

figure(2);
title('Coarse');
plot(real(rx_coarse));

%%% Perform fine synchronization

[rx_fine,phErr] = fine(rx_coarse);

figure(3);
title('Fine');
plot(real(rx_fine));

end
