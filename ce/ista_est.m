function [H_ista, Info] = ista_est(H_noisy, GenPar, ChanInfo)
% ista based channel estimation
    Info = [];    
    H_ls_comb = H_noisy(1:GenPar.comb:end,:);
    
    comb = GenPar.comb;
    Nsc_comb = GenPar.Nsc / GenPar.comb;
    Nfft = ChanInfo.Nfft;

    %Nfft_us = Nfft + 16;
    num_iter = GenPar.ista_iter;
    
    A_dft = dftmtx(Nfft) ./ sqrt(Nfft);
    %A_dft = dftmtx(Nfft) ./ sqrt(Nsc_comb);
    A_dft_int = dftmtx(Nfft*GenPar.comb) ./ sqrt(Nfft);
    %A_dft_int = dftmtx(Nfft*GenPar.comb) ./ sqrt(Nsc_comb);
    
    % init determine window sizes
    win_max = fix( ChanInfo.tau_max * GenPar.delta_f * Nfft * GenPar.comb );
    win_min = fix( ChanInfo.win_guard * Nfft * GenPar.comb / 512);
    
    pos_signal = [(1:win_max), (Nfft - win_min + 1:Nfft)];
    pos_sig_int = pos_signal;
    
    pos_noise = [win_max+1:Nfft-win_min];
    
    pos_high = find(pos_sig_int > Nfft/2);
    pos_sig_int(pos_high) = pos_sig_int(pos_high) - Nfft;
    x_int = (0:GenPar.Nsc-1).';
    x_int = x_int ./ comb;
    A_sig_int = exp( -1j * 2*pi * x_int * (pos_sig_int-1) ./ Nfft ) ./ sqrt(Nfft);
    
    A_sig = A_dft(1:Nsc_comb, pos_signal);
    
    A_noise = A_dft(1:Nsc_comb, pos_noise);
    win_all = length(pos_signal);
    
    R = H_ls_comb;
    X_ista = zeros(win_all, GenPar.Nrx);
    
    %sigma_noise = GenPar.sigma_noise_sc;
    
    sigma_est_noise = A_noise'*H_ls_comb;
    sigma_est_noise = mean( mean(abs(sigma_est_noise).^2, 2) );
    
    sigma_noise = sigma_est_noise;
    
    
    %nu_arr = [1.0:-(1.0/num_iter):0.0];
    nu_arr = ones(1, num_iter);
    
    if ~GenPar.do_amp
        for iter_idx = 1:num_iter
            CorrMx = A_sig'*R;
            X_ista = X_ista + CorrMx;

            lambda = nu_arr(iter_idx) * sigma_noise;
            [X_filt, num_non_zeros] = soft_thr_mmv(X_ista, lambda);
            X_ista = X_filt;

            R = H_ls_comb - A_sig * X_filt;
            %norm(R)
            fprintf('Residual at iter=%d equal =%f\n', iter_idx, norm(R));
        end
    else
        Z = H_ls_comb;
        N_dict = size(A_sig, 1);
        num_non_zeros_priv = 0.0;
        X_priv = zeros(win_all, GenPar.Nrx);
        
        iter_warm = GenPar.ista_iter_warm;
        X_iter = {};
        
        for iter_idx = 1:num_iter
            CorrMx = A_sig'*Z;
            X_ista = X_ista + CorrMx;

            if iter_idx > iter_warm
                lambda = ( norm(X_ista(:) - X_iter{iter_idx-1}(:))^2 ./ length(X_ista(:)) );
            else
                lambda = nu_arr(iter_idx) * sigma_noise;
            end
            [X_filt, num_non_zeros] = soft_thr_mmv(X_ista, lambda);
            X_ista = X_filt;

            X_iter{iter_idx} = X_ista;
            
            Z = H_ls_comb - A_sig * X_ista + 0.5*(num_non_zeros/N_dict) * Z;
            X_priv = X_ista;
            num_non_zeros_priv = num_non_zeros;
            %norm(R)
            fprintf('Residual at iter=%d equal =%f\n', iter_idx, norm(Z));
        end
    end
    
    H_ista = A_sig * X_filt;
    H_ista = A_sig_int * X_filt;
    
    %err1 = norm(A_sig - A_sig_int)
end

function [X_filt, num_non_zeros] = soft_thr_mmv(X_ista, lambda)
    
    [WinAll, Nrx] = size(X_ista);   
    pdp = mean(abs(X_ista).^2, 2);
    pos_high = find(pdp > lambda);
    num_non_zeros = length(pos_high);
    filterIsta = zeros(WinAll, 1);
    filterIsta(pos_high) = (pdp(pos_high) - lambda) ./ pdp(pos_high);
    FiltArr = repmat(filterIsta, 1, Nrx);
    X_filt = X_ista .* FiltArr;

end
