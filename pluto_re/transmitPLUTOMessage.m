%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Author   : David Marcus, PE
%%% Date     : 03 NOV 2019
%%% Rev      : 1 - Initial Release
%
%%%%%%%%%%%%%%% Transmit PLUTO Message %%%%%%%%%%%%%%%
%
%%%%% Inputs %%%%%
%
%%% If inputs are empty, then defaults will be used. 
%
%%% center_frequency:   RF Frequency transmit frequency (default = 915Mhz)
%%% sample_rate:        Data rate in symbols/second times oversamplerate
%%%                     in samples/symbol (default = 16 samples/second)
%%% data_rate:          Bit rate in samples/second
%%% message_num:        Message number to transmit
%
%%%%% Outputs %%%%%
%
%%% plutoTX:            ADALM PLTUO system object
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [plutoTx, message] = transmitPLUTOMessage(PrCfg, message_num, center_frequency,sample_rate,...
    data_rate)

global plutoRadio

if ~exist('message_num','var')
    message_num = 1; % Set to 1, 2, or 3
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

[IQmessage ,message_length, message] = createMessage(message_num, PrCfg);

sps = sample_rate/data_rate;

txFilter = comm.RaisedCosineTransmitFilter( ...
    'OutputSamplesPerSymbol',sps,'FilterSpanInSymbols',2);

tx_signal = txFilter(IQmessage);

plutoRadio = findPlutoRadio;

plutoTx = plutoTxInit(center_frequency,sample_rate);

transmitRepeat(plutoTx,tx_signal);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
