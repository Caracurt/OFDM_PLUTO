%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Author   : David Marcus, PE
%%% Date     : 03 NOV 2019
%%% Rev      : 1 - Initial Release  
%
%%%%%%%%%%%%%%% RTL-SDR Receiver Initialization %%%%%%%%%%%%%%%
%
%%%%% Inputs %%%%%
%
%%% If inputs are empty, then defaults will be used. 
%
%%% center_frequency:   RF Frequency transmit frequency (default = 915Mhz)
%%% sample_rate:        Data rate in symbols/second times oversamplerate
%%%                     in samples/symbol (default = 300e3 samples/second)  
%%% samples_per_frame:  Total number of samples to capture
%
%%%%% Outputs %%%%%
%
%%% RTLSDRrx             RTL-SDR system object
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RTLSDRrx = RTLSDRRxInit(center_frequency,sample_rate,...
    samples_per_frame)

%%% Set input variables to default if not input with function
if ~exist('center_frequency','var')
    center_frequency = 915e6;
end
if ~exist('sample_rate','var')
    sample_rate = 300e3;
end

%%% RTLSDR Radio Parameters %%%

RTLSDRrx = comm.SDRRTLReceiver('0','CenterFrequency',center_frequency,...
    'SampleRate',sample_rate,'SamplesPerFrame',samples_per_frame,...
    'OutputDataType','double');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%