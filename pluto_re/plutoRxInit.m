%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Author   : David Marcus, PE
%%% Date     : 03 NOV 2019
%%% Rev      : 1 - Initial Release  
%
%%%%%%%%%%%%%%% Pluto Receiver Initialization %%%%%%%%%%%%%%%
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
%%% plutoRx             ADALM PLUTO system object
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plutoRx = plutoRxInit(center_frequency,sample_rate,...
    samples_per_frame)

%global plutoRadio

%%% Set input variables to default if not input with function
if ~exist('center_frequency','var')
    center_frequency = 915e6;
end
if ~exist('sample_rate','var')
    sample_rate = 300e3;
end
                                                               
name = 'Pluto';
%plutoRadio = findPlutoRadio;
%%% Pluto Receiver Radio Parameters %%%

plutoRx = sdrrx(name,'RadioID','ip:192.168.2.2',...
    'BasebandSampleRate',sample_rate,'CenterFrequency',center_frequency,...
    'OutputDataType','double','SamplesPerFrame',samples_per_frame);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%