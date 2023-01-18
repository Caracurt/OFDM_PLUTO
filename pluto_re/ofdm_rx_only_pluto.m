% transmit and receive harmonic function on sdr
if exist('plutoTx')
    release(plutoTx)
end

if exist('plutoRx')
    release(plutoRx)
end

close all; clear all; clc;

%global plutoRadio

center_frequency = 2.4e9;
sample_rate = 8 * 940e3;

%data_freq = sample_rate/10; % frequecny of sin
%data_rate = 

% data_length = 128;
% guard_length = 128;
% freq_data = 0/data_length;
% x_arr = 0:data_length-1;
% x_arr = x_arr';
% sig_tx = exp(1j*2*pi*x_arr*freq_data);
% guard_int = zeros(guard_length, 1);
% %tx_mess = complex( [sig_tx; guard_int; -sig_tx; guard_int] );
% tx_mess = complex( [sig_tx; guard_int; -sig_tx; guard_int; guard_int; guard_int] );
% 
% %tx_mess = complex( [sig_tx; sig_tx; sig_tx; sig_tx] );
% 
% sig_len = length(tx_mess);
% 
% figure(10);
% plot(real(tx_mess))

% Pluto transmit

% prepare message 

rng(1);

% OFDM message
delta_t = 1 / sample_rate;

delta_f = 7.5e3; % subcarrier spacing
frac_guard = 0.5;

N_sc_av = fix( sample_rate / delta_f );
N_fft = 2^( fix(log2(N_sc_av)) );
N_sc_use = fix(N_fft * frac_guard);
guard_length = fix( 0.2 * N_fft );
%guard_length = 0;
CP_len = fix( N_fft * 0.2 );
num_bits_sym = 1;

um_bits_sym = 1;
BlockSize = N_sc_use * num_bits_sym; 

N_repeat = 10; % number of repeats for tx

% create preambule
bark = comm.BarkerCode('SamplesPerFrame',N_fft/2,'Length',13);
preamble = complex(bark());
release(bark);

preamble_full = [preamble; preamble];
preamble_full_cp = [preamble_full(end - CP_len + 1:end); preamble_full];

% pilot frame
data_stream = randi([0, 1], BlockSize, 1);
mod_sym_pilot = 1.0 - 2.0 .* data_stream; % BPSK modulation

tx_ofdm_sym = zeros(N_fft, 1);
tx_ofdm_sym(2:N_sc_use/2 + 1) = mod_sym_pilot(N_sc_use/2 + 1 : end);
tx_ofdm_sym(end - N_sc_use/2 + 1 : end) = mod_sym_pilot(1:N_sc_use/2);

time_ofdm_sym_pilot = ifft(tx_ofdm_sym);
time_ofdm_sym_cp_pilot = [time_ofdm_sym_pilot(end - CP_len + 1:end); time_ofdm_sym_pilot];


% data frame
data_stream = randi([0, 1], BlockSize, 1);
mod_sym = 1.0 - 2.0 .* data_stream; % BPSK modulation

tx_ofdm_sym = zeros(N_fft, 1);
tx_ofdm_sym(2:N_sc_use/2 + 1) = mod_sym(N_sc_use/2 + 1 : end);
tx_ofdm_sym(end - N_sc_use/2 + 1 : end) = mod_sym(1:N_sc_use/2);

time_ofdm_sym = ifft(tx_ofdm_sym);
time_ofdm_sym_cp = [time_ofdm_sym(end - CP_len + 1:end); time_ofdm_sym];

guard_frame = zeros(guard_length, 1);

full_frame = [preamble_full_cp; time_ofdm_sym_cp_pilot; time_ofdm_sym_cp; guard_frame];
FrameSize = length(full_frame);
sig_len = length(full_frame);

% upsampling and RRC
%sps = sample_rate/data_rate;

%txFilter = comm.RaisedCosineTransmitFilter( ...
%    'OutputSamplesPerSymbol',sps,'FilterSpanInSymbols',2);

%tx_signal = txFilter(IQmessage);

%plutoTx = plutoTxInit(center_frequency,sample_rate);
%transmitRepeat(plutoTx,full_frame);
% actual Tx

% pause
%pause(30)

%return;

samples_per_frame = 20*sig_len;

%RTLSDRrx = RTLSDRRxInit(center_frequency,sample_rate,samples_per_frame);
%rx_signal = RTLSDRrx();
%release(RTLSDRrx);

%plutoRx = plutoRxInit(center_frequency,sample_rate,samples_per_frame);
%rx_signal = plutoRx();
%release(plutoRx);

load('rx_sig_m1.mat');
%rx_signal = rx_sig(:);
rx_signal = rx_signal(1:38880);

figure(20);
plot(10*log10(abs(rx_signal)));


%frame_rx = frame_rx_cfo(init_pos + 1: init_pos + reception_win);
frame_rx = rx_signal;
% step 1 - search begging of the frame
len_rec = length(frame_rx);
corr_len = N_fft / 2;

metric_p = [];

for start_pos = 0:len_rec-N_fft-1
    
    A1 = frame_rx(start_pos + 1: start_pos + N_fft/2);
    A2 = frame_rx(start_pos + N_fft/2 + 1: start_pos + N_fft);
    
    corr = A1'*A2;
    metric_p = [metric_p; corr];
    
end

filt_1 = ones(CP_len, 1);
metric_p1 = conv(metric_p, filt_1, 'same');


figure(1);
plot(abs(metric_p));
grid on

figure(2);
plot(abs(metric_p1));
grid on

[max_val, pos_start] = max(abs(metric_p1(1:FrameSize)));
pos_start = pos_start + fix(CP_len / 2);
%pos_start = CP_len / 2;

cfo_est = angle( metric_p(pos_start) ) / (2*pi) / (N_fft/2);
%delta_cfo = cfo_idl - cfo_est;

x_arr = (0:length(frame_rx)-1).';
%cfo_est = 0.0;

for cf_flag = [0, 1]

    if (cf_flag == 0)
        cfo_est_use = 0;
    else
        cfo_est_use = cfo_est;
    end
    cfo_corr = exp(1j * (-2*pi*cfo_est_use) * x_arr);
    frame_rx_cfo = frame_rx .* cfo_corr;


    frame_ta = frame_rx_cfo( pos_start : pos_start + N_fft + (CP_len + N_fft) + (CP_len + N_fft) );

    % compensate CFO
    %x_arr = (0:length(frame_ta)-1).';
    %cfo_corr = exp(1j * (-cfo_est) * x_arr);
    %frame_ta = frame_ta .* cfo_corr;

    pilot_frame = frame_ta(N_fft + CP_len + 1 : N_fft + CP_len + N_fft);
    pilot_freq = fft(pilot_frame);


    pilot_freq_bb = [pilot_freq(end - N_sc_use/2 + 1 : end) ; pilot_freq(2 : N_sc_use/2 + 1)];
    %mod_sym_pilot

    h_ls = pilot_freq_bb ./ mod_sym_pilot;
    
    % Channel estimation
    win_for = 70;
    win_back = 50;
    
    N_sc_ls = length(h_ls);
    N_fft_ce = 2^( ceil(log2(N_sc_ls) ) + 2 );
    h_ls_zp = zeros(N_fft_ce, 1);
    h_ls_zp(1:N_sc_ls) = h_ls;
    
    % rebuilding
    reb_right = 6;
    reb_left = 6;
    
    do_reb = true;
    
    
    if do_reb
        % right
        h_ls_zp(N_sc_ls + 1 : N_sc_ls + reb_right) = h_ls_zp(N_sc_ls);

        %left
        h_ls_zp(end - reb_left + 1 : end) = h_ls_zp(1);
    end
    
    h_time = ifft(h_ls_zp) * sqrt(N_fft_ce);
    
    % noise est;
    h_noise = h_time(win_for+1:end-win_back);
    sigma_noise = mean(abs(h_noise).^2);
    
    filt_hw = zeros(N_fft_ce, 1);
    filt_hw(1:win_for) = 1.0;
    filt_hw(N_fft_ce - win_back + 1 : end) = 1.0;
    
    h_time_filt = h_time .* filt_hw;
    
    % soft win
    pdp_time = abs(h_time_filt).^2;
    pos_high = find(pdp_time  > sigma_noise);
    filt_sw = zeros(N_fft_ce, 1);
    filt_sw(pos_high) = ( pdp_time(pos_high) - sigma_noise) ./ pdp_time(pos_high);
    
    % sw
    h_time_filt = h_time_filt .* filt_sw;
    
    h_freq_us = fft(h_time_filt) ./ sqrt(N_fft_ce);    
    h_ls = h_freq_us(1:N_sc_ls);
    % chane est end
    
    figure(3 + cf_flag*10);
    plot(abs(pilot_freq));

    figure(4 + cf_flag*10);
    plot(real(h_ls));
    title('LS channel estimate');
    grid on;
    
    figure(5 + cf_flag*10);
    plot(abs(ifft(h_ls)));
    title('LS channel estimate: time');
    grid on;

    % data frame
    data_frame = frame_ta(N_fft + CP_len + N_fft + CP_len + 1: N_fft + CP_len + N_fft + CP_len + N_fft);
    data_freq = fft(data_frame);


    data_freq_bb = [data_freq(end - N_sc_use/2 + 1 : end) ; data_freq(2 : N_sc_use/2 + 1)];

    %mod_sym

    demod_sym = data_freq_bb ./ h_ls;

    %demod_sym = [demod_sym(N_sc_use/2 + 1 : end); demod_sym(1:N_sc_use/2)];

    evm_est = norm(demod_sym - mod_sym) ./ norm(mod_sym);

    figure(6 + cf_flag*10);
    scatterplot(demod_sym);
    grid on;
    title(['CFOcomp=', num2str(cf_flag)])


    bit_dec = double(demod_sym < 0);
    ber = sum(bit_dec ~= data_stream) / length(data_stream);

    fprintf('BER=%f CFOcomp=%d EVMest=%f\n', ber, cf_flag , evm_est)
end


if 0
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
end

