function [H_sw, Info] = sw_est(H_noisy, GenPar, ChanInfo)
% sw channel estimation
    Info = [];    
    H_ls_comb = H_noisy(1:GenPar.comb:end,:);
    
    Nsc_comb = GenPar.Nsc / GenPar.comb;
    
    if 0
        if GenPar.do_pow2
            Nfft = 2^ceil(log2(Nsc_comb));
        else
            Nfft = Nsc_comb + 16;
        end
    end
    Nfft = ChanInfo.Nfft;
    
    A_dft = dftmtx(Nfft) ./ sqrt(Nfft);
    A_dft_int = dftmtx(Nfft*GenPar.comb) ./ sqrt(Nfft);
    
    % init determine window sizes
    %win_max = fix( ChanInfo.tau_max * GenPar.delta_f * Nfft * GenPar.comb  );
    win_max = ChanInfo.win_max;
    %win_min = fix( ChanInfo.win_guard * Nfft * GenPar.comb / 512);
    win_min = ChanInfo.win_min;
    
    pos_signal = [(1:win_max), (Nfft - win_min + 1:Nfft)];
    W_sig = length(pos_signal);
    pos_noise = (win_max+1:Nfft-win_min);
    W_noise = length(pos_noise);
    
    % dictionary extraction
    A_sig = A_dft(:, pos_signal);
    A_noise = A_dft(:, pos_noise);
    
    % upsampling
    H_ls_us = zeros(Nfft, GenPar.Nrx);
    H_ls_us(1:Nsc_comb, :) = H_ls_comb;
    
    % let us make rebuilding
    if GenPar.do_rebuild
        [H_left, H_right] = calc_rebuild(H_ls_comb, GenPar, ChanInfo);

        H_ls_us(Nsc_comb+1:Nsc_comb+GenPar.reb_num,:) = H_right;
        H_ls_us(Nfft:-1:Nfft-GenPar.reb_num+1,:) = H_left;
    end
    
    
    
    if 0
    ht = H_ls_us(:,1);
    ht = circshift(ht, [6, 0]);
    figure(42);
    plot(real(ht), '-r');
    keyboard;
    end
    
    H_time = A_sig'*H_ls_us;
    H_noise = A_noise'*H_ls_us;
    
    pdp = mean(abs(H_time).^2, 2);
    sigma_noise = mean( mean(abs(H_noise).^2 , 2) );
    
    pos_high = find(pdp > sigma_noise);
    filtPdp = zeros(W_sig, 1);
    
    filtPdp(pos_high) = (pdp(pos_high) - sigma_noise) ./ (pdp(pos_high));
    
    FiltPdp = repmat(filtPdp, 1, GenPar.Nrx);
    H_time_filt = H_time .* FiltPdp;
    
    ht = zeros(Nfft, GenPar.Nrx);
    ht(pos_signal,:) = H_time_filt;
    
    % interpolation
    H_time_us = zeros(Nfft*GenPar.comb, GenPar.Nrx);
    H_time_us(1:Nfft/2,:) = ht(1:Nfft/2,:);
    H_time_us(end - Nfft/2 + 1:end,:) = ht(Nfft/2 + 1:end,:);
    
    H_freq_us = A_dft_int * H_time_us;
    
    H_sw = H_freq_us(1:GenPar.Nsc,:);
end

