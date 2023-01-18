%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Author   : David Marcus, PE
%%% Date     : 03 NOV 2019
%%% Rev      : 1 - Initial Release
%
%%%%%%%%%%%%%%% Receive RTL-SDR Message %%%%%%%%%%%%%%%
%
%%%%% Inputs %%%%%
%
%%% If inputs are empty, then defaults will be used. 
%
%%% center_frequency:   RF Frequency transmit frequency (default = 915Mhz)
%%% sample_rate:        Data rate in symbols/second times oversamplerate
%%%                     in samples/symbol (default = 16 samples/second)
%%% data_rate:          Bit rate in samples/second
%%% message_num:        Message number to decode      
%
%%%%% Outputs %%%%%
%
%%% decodedMessage:     Decoded message from receiver
%%% rx_signal:          Samples received from RTL-SDR radio
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [decodedMessage,rx_signal] = receiveRTLSDRMessage(PrCfg, message_num, center_frequency,...
    sample_rate,data_rate)

if ~exist('message_num','var')
    message_num = 1;
end

if ~exist('center_frequency','var')
    center_frequency = 915e6;
end

if ~exist('sample_rate','var')
    sample_rate = 300e3;
end

if ~exist('data_rate','var')
    data_rate = 20e3;
end

if message_num == 1
    message = 'www.EEmaginations.com';
elseif message_num == 2
    message = 'David E. Marcus, PE';
elseif message_num == 3
    message = 'Hello Chunia !!!';
elseif message_num == 4
    message = 'Yana, sweet girl!';
end

message_length = length(message);
guard_length = PrCfg.guard_len;
preamble_length = PrCfg.preamb_len * 2;

sps = sample_rate/data_rate;

%%% Capture about 10 frames of data

samples_per_frame = 10*sps*(guard_length + message_length*8 + preamble_length);

RTLSDRrx = RTLSDRRxInit(center_frequency,sample_rate,samples_per_frame);

rx_signal = RTLSDRrx();

figure(1);
plot(10*log10(abs(rx_signal)))
figure(2);
plot((real(rx_signal)))
release(RTLSDRrx);

rxFilter = comm.RaisedCosineReceiveFilter(...
    'InputSamplesPerSymbol',sps,'DecimationFactor',sps,...
    'FilterSpanInSymbols',2);

decodedMessage = decodeMessage(PrCfg, rxFilter(rx_signal),message_length,sample_rate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%