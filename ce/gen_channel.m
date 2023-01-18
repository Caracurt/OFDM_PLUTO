function [H_out, Info] = gen_channel(GenPar, ChanInfo);
% generate channel
    % input structs
% GenPar.Nsc = 192;
% GenPar.comb = 2;
% GenPar.delta_f = 30e3; % supcarrier spacing
% GenPar.Nrx = 16;
% 
% % channel parameters
% ChanInfo.num_taps = 8;
% ChanInfo.tau_max = 5e-6;
% 
% ChanInfo.win_max = fix( tau_max * delta_f * Nsc  );
% ChanInfo.win_guard = 24; % max sync error for Nfft=2048
% ChanInfo.win_min = fix( win_guard * Nsc / 2048);
    
    % generate tap frequencies
    tau_rms = ChanInfo.tau_max/4;
    tau_abs = abs( tau_rms .* randn(ChanInfo.num_taps, 1));
    
    % generate amplitudes
    A_arr = zeros(ChanInfo.num_taps, GenPar.Nrx);
    scaling_mask = repmat( exp(-tau_abs./tau_rms),1 , GenPar.Nrx );
    
    A_arr = (1.0/sqrt(2)) * ( randn(ChanInfo.num_taps, GenPar.Nrx) + 1j*randn(ChanInfo.num_taps, GenPar.Nrx) );
    A_arr = scaling_mask .* A_arr;
    
    x_freq = (0:1:GenPar.Nsc-1).*GenPar.delta_f;
    Exp_Arr = exp(-1j*2*pi*(x_freq(:)*tau_abs.'));
    
    Info.freq_out = tau_abs * GenPar.delta_f;
    
    H_out = Exp_Arr * A_arr;
    
    H_out = H_out ./ norm(H_out(:)) .* sqrt(GenPar.Nrx);
    %Info = [];
end

