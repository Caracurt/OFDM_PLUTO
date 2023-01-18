%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Author   : David Marcus, PE
%%% Date     : 03 NOV 2019
%%% Rev      : 1 - Initial Release
%
%%%%%%%%%%%%%%% Create Message %%%%%%%%%%%%%%%
%
%%%%% Inputs %%%%%
%
%%% message_num:    message_num = 1 = www.EEmaginations.com (default)
%%%                 message_num = 2 = David E. Marcus, PE
%%%                 message_num = 3 = Hello World !!!
%
%%%%% Outputs %%%%%
%
%%% IQMessage:  BPSK modulated, frame structured as follows:
%
%%% Frame = |Guard|Preamble|Message|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [IQMessage,message_length, message] = createMessage(message_num, PrCfg)

if (~exist('message_num','var')) || (message_num == 1)
    message = 'www.EEmaginations.com';
elseif message_num == 2
    message = 'David E. Marcus, PE';
elseif message_num == 3
    message = 'Hello Chunia !!!';
elseif message_num == 4
    message = 'Yana, sweet girl!';
end

message_length = length(message);

%%% Generate Barker Code for Preamble

bark = comm.BarkerCode('SamplesPerFrame',PrCfg.preamb_len,'Length',13);
preamble = complex(bark());
release(bark);

%%% Convert message to bits (1s and 0s)
index = 1;
message_binary = zeros(length(message)*8,1);
for i = 1:8:length(message_binary)-7
    message_bits_string = dec2bin(message(index),8)';
    for j = 1:8
        message_binary(i+j-1) = str2double(message_bits_string(j));
    end
    index = index + 1;
end

%%% Modulate data to BPSK (Symbol Mapping: 1 -> 1, 0 -> -1)

BPSKmodulator = comm.BPSKModulator('PhaseOffset',pi);

BPSKdata = BPSKmodulator(message_binary);
release(BPSKmodulator);

%guard = complex(0.0001*randn(1,PrCfg.guard_len)');
guard = complex(0.000*randn(1,PrCfg.guard_len)');

IQMessage = [guard;preamble;preamble;BPSKdata];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%