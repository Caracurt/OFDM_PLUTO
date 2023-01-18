function [H_freq, Info] = omp_est(H_noisy, GenPar, ChanInfo)
    % OMP channel estimate
    Info = [];    
    H_ls_comb = H_noisy(1:GenPar.comb:end,:);
    comb = GenPar.comb;
    
    Nsc_comb = GenPar.Nsc / GenPar.comb;
    Nfft = ChanInfo.Nfft;

    %Nfft_us = Nfft + 16;
    num_iter = GenPar.num_omp_iter;
    A_dft = dftmtx(Nfft) ./ sqrt(Nfft);
    %A_dft_int = dftmtx(Nfft*GenPar.comb*factor_us) ./ sqrt(Nfft*factor_us);
    
    % init determine window sizes
    win_max = fix( ChanInfo.tau_max * GenPar.delta_f * Nfft * GenPar.comb );
    win_min = fix( ChanInfo.win_guard * Nfft * GenPar.comb / 512);
    
    pos_signal = [(1:win_max), (Nfft - win_min + 1:Nfft)];
        
    A_sig = A_dft(1:Nsc_comb, pos_signal);
    
    [H_time, pos, H_denoise, A_pos] = omp_core(H_ls_comb, A_sig, num_iter, GenPar.sigma_noise_sc);
    
    pos_sig_int = pos_signal(pos);
    pos_high = find(pos_sig_int > Nfft/2);
    pos_sig_int(pos_high) = pos_sig_int(pos_high) - Nfft;
    x_int = (0:GenPar.Nsc-1).';
    x_int = x_int ./ comb;
    A_sig_int = exp( -1j * 2*pi * x_int * (pos_sig_int-1) ./ Nfft ) ./ sqrt(Nfft);
    
    %norm(A_pos - A_sig_int)
    H_freq = A_sig_int * H_time;
    
    if 0
        % interpolation
        A_1 = dftmtx(Nfft) ./ sqrt(Nfft);
        A_2 = dftmtx(Nfft*GenPar.comb) ./ sqrt(Nfft);
        H_time = A_1'*H_denoise;

        H_us = [H_time(1:Nfft/2,:); zeros( (GenPar.comb-1)*Nfft, GenPar.Nrx ); H_time(Nfft/2 + 1:end,:)];

        H_freq = A_2 * H_us;
    end
end

function [H_time_tmp, pos, H_denoise, A_cut] = omp_core(H_ls_comb, A_sig, num_iter, sigma_noise_sc)

    [N_sc, N_rx] = size(H_ls_comb);
    [~, N_time] = size(A_sig);

    H_time = zeros(N_time, N_rx);
    R = H_ls_comb;
    
    pos = [];
    for i = 1:num_iter
        r = A_sig'*R;
        r = mean( abs(r).^2, 2);
        [~, pos_max] = max(r);
        
        pos = [pos, pos_max];
        
        A_cut = A_sig(:, pos);
        H_time_tmp = pinv(A_cut'*A_cut) * A_cut' * H_ls_comb;
        H_denoise = A_cut * H_time_tmp;
        
        R = H_ls_comb - H_denoise;
        if (norm(R(:))^2 ./ length(R(:)) < 1.0*sigma_noise_sc)
            break;
        end
        rr = norm(H_ls_comb - H_denoise);
    end

    %pdp_x = mean(abs(H_time_tmp).^2, 2);
    %H_time_tmp = pinv(A_cut'*A_cut + diag(sigma_noise_sc./pdp_x)) * A_cut' * H_ls_comb;
    H_denoise = A_cut * H_time_tmp;
    
end