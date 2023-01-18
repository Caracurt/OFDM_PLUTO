%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Author   : David Marcus, PE
%%% Date     : 03 NOV 2019
%%% Rev      : 1 - Initial Release
%
%%%%%%%%%%%%%%% Decode Message %%%%%%%%%%%%%%%
%
%%%%% Inputs %%%%%
%
%%% rxSignal:       Received RF signal
%%% message_length: number of characters in message
%%% fs:             Sample rate (Hz)
%
%%%%% Outputs %%%%%
%
%%% decodedMessage: Decoded message from receiver
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [decodedMessage,rx_coarse,rx_fine,data_start_index] = ...
    decodeMessage(PrCfg, rxSignal,message_length,fs)

if ~exist('fs','var')
    fs = 300000;
end

%%% Generate Barker Code for Preamble

bark = comm.BarkerCode('SamplesPerFrame',PrCfg.preamb_len,'Length',13);
preamble = bark();
release(bark);

preamble_long = [preamble;preamble];

%%% Set up objects to correct received signal

coarse = comm.CoarseFrequencyCompensator('SampleRate',fs, ...
    'Modulation','BPSK','FrequencyResolution',1);
fine = comm.CarrierSynchronizer( ...
    'SamplesPerSymbol',1,'Modulation','BPSK');
prbdet = comm.PreambleDetector([preamble;preamble],'Threshold',3,'Detections','All');

%%% Perform coarse synchronization

[rx_coarse,estfreqoff] = coarse(rxSignal);

%%% Perform fine synchronization

[rx_fine,phErr] = fine(rx_coarse);

%%% Find data start using preamble

[idx1,detmet1] = prbdet(rx_fine);

sorted_detmet1 = sort(detmet1,'descend');

prbdet.Threshold = floor(sorted_detmet1(2));

[idx2,detmet2] = prbdet(rx_fine);

data_start_index = idx2(1) + 1;
data_end_index = data_start_index + message_length*8 - 1;

%%% guard length = 100
guard_start_index = idx2(1) - PrCfg.guard_len - length(preamble_long) + 1;
guard_end_index = guard_start_index + PrCfg.guard_len - 1;

preamble_start_index = idx2(1) - length(preamble_long) + 1;
preamble_end_index = preamble_start_index + length(preamble_long) - 1;

%%% Check if synchronization algorithm added 180 degree phase offset

rx_preamble_symbols = rx_fine(preamble_start_index:preamble_end_index);
rx_preamble_bits = real(rx_preamble_symbols);
rx_preamble_data = zeros(length(preamble_long),1);
for i = 1:length(preamble_long)
    if rx_preamble_bits(i) > 0
        rx_preamble_data(i) = 1;
    else
        rx_preamble_data(i) = -1;
    end
end

phOffsetEst = angle(preamble_long.*rx_preamble_data);
ph_offset_degrees = rad2deg(mean(phOffsetEst));

if abs(ph_offset_degrees) > 90
    rx_fine = -rx_fine;
end

%%% Begin decoding algorithm 

rx_message_symbols = rx_fine(data_start_index:data_end_index);
rx_message_bits = real(rx_message_symbols);
rx_message_data = rx_message_bits > 0;

rx_message_data_string = num2str(rx_message_data);
rx_message_characters = strings(message_length,1);
index = 1;
for i = 1:8:length(rx_message_data)-7
    rx_message_characters(index) = strcat(rx_message_data_string(i:i+7))';
    index = index + 1;
end

decodedMessage = char(bin2dec(rx_message_characters))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%