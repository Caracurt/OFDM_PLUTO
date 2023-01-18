
center_frequency = 1.2e9;
sample_rate = 1000e3;
samples_per_frame = 4096;
RTLSDRrx = RTLSDRRxInit(center_frequency,sample_rate,samples_per_frame);
rx_signal = RTLSDRrx();

figure(1);
plot(10*log10(abs(fft(rx_signal))))
%figure(2);
%plot((real(rx_signal)))

release(RTLSDRrx);