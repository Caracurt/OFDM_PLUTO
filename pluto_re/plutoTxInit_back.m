%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Author   : David Marcus, PE
%%% Date     : 03 NOV 2019
%%% Rev      : 1 - Initial Release  
%
%%%%%%%%%%%%%%% Pluto Transmitter Initialization %%%%%%%%%%%%%%%
%
%%%%% Inputs %%%%%
%
%%% If inputs are empty, then defaults will be used. 
%
%%% center_frequency:   RF Frequency transmit frequency (default = 915Mhz)
%%% sample_rate:        Data rate in symbols/second times oversamplerate
%%%                     in samples/symbol (default = 300e3 samples/second)        
%
%%%%% Outputs %%%%%
%
%%% plutoTx:            ADALM PLUTO system object
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plutoTx = plutoTxInit(center_frequency,sample_rate)

global plutoRadio

%%% Set input variables to default if not input with function
if ~exist('center_frequency','var')
    center_frequency = 915e6;
end
if ~exist('sample_rate','var')
    sample_rate = 300e3;
end

%%% Set Pluto radio parameters

% configurePlutoRadio('AD9364');   % Configure as AD9364 if need larger RF
                                 % RF Bandwidth = 70MHz to 6GHz
                                 
name = 'Pluto';
gain = 0;

%%% Pluto Radio Parameters %%%
plutoTx = sdrtx(name,'RadioID',plutoRadio.RadioID,'CenterFrequency',...
    center_frequency,'BasebandSampleRate',sample_rate,'Gain',gain);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%