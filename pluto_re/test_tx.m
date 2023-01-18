%addpath('./SDR Functions');

if exist('plutoTx')
    release(plutoTx)
end

close all; clear all; clc;

PrCfg = struct;
PrCfg.guard_len = 500;
PrCfg.preamb_len = 2^7;

center_frequency = 1300e6; % 915e6 1.2e9
%sample_rate = 940e3; % 300e3
sample_rate = 940 * 1e3; % 300e3
data_rate = (sample_rate/10); % 20e3
mess_idx = 4;
[plutoTx, message] = transmitPLUTOMessage(PrCfg, mess_idx, ...
    center_frequency, sample_rate, data_rate);

pause(0.1)

% rtl sdr
%[decodedMessage,rx_signal] = receiveRTLSDRMessage(PrCfg, mess_idx, ...
%    center_frequency, sample_rate, data_rate);

% adlm pluto
[decodedMessage,rx_signal] = receivePLUTOMessage(PrCfg, mess_idx, ...
   center_frequency, sample_rate, data_rate);

release(plutoTx)

[decodedMessage , " ", message]