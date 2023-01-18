function [H_imp, Info] = timp_est(H_noisy, GenPar, ChanInfo)
% imp channel estimation
    Info = [];    
    H_ls_comb = H_noisy(1:GenPar.comb:end,:);
    
    comb = GenPar.comb;
    Nsc_comb = GenPar.Nsc / GenPar.comb;
    Nfft = ChanInfo.Nfft;

    num_iter = GenPar.ista_iter;
    
    A_dft = dftmtx(Nfft) ./ sqrt(Nsc_comb);
    A_dft_int = dftmtx(Nfft*GenPar.comb) ./ sqrt(Nsc_comb);
    
    % init determine window sizes
    win_max = fix( ChanInfo.tau_max * GenPar.delta_f * Nfft * GenPar.comb ) + 2;
    win_min = fix( ChanInfo.win_guard * Nfft * GenPar.comb / 512);
    
    pos_signal = [(1:win_max), (Nfft - win_min + 1:Nfft)];
    pos_sig_int = pos_signal;
    
    pos_high = find(pos_sig_int > Nfft/2);
    pos_sig_int(pos_high) = pos_sig_int(pos_high) - Nfft;
    x_int = (0:GenPar.Nsc-1).';
    x_int = x_int ./ comb;
    A_sig_int = exp( -1j * 2*pi * x_int * (pos_sig_int-1) ./ Nfft ) ./ sqrt(Nsc_comb);
    
    A_sig = A_dft(1:Nsc_comb, pos_signal);   
    
    A_sig_us = A_dft(:, pos_signal);
    
    win_all = length(pos_signal);
    
    pos_noise = [win_max+1:Nfft-win_min];
    A_noise = A_dft(1:Nsc_comb, pos_noise);
    
    sigma_est_noise = A_noise'*H_ls_comb;
    sigma_est_noise = mean( mean(abs(sigma_est_noise).^2, 2) );
    
    sigma_noise = sigma_est_noise;
    
    S = A_sig'*A_sig;
    S = S * diag(  1.0 ./ diag(S) ) * 1.0;
    
    %S = A_sig_us'*A_sig_us;
    %S = S * diag(  1.0 ./ diag(S) ) * 1.0;
    
    % let us make rebuilding
    H_ls_us = zeros(Nfft, GenPar.Nrx);
    H_ls_us(1:Nsc_comb, :) = H_ls_comb;
    
    if GenPar.do_rebuild
        [H_left, H_right] = calc_rebuild(H_ls_comb, GenPar, ChanInfo);

        H_ls_us(Nsc_comb+1:Nsc_comb+GenPar.reb_num,:) = H_right;
        H_ls_us(Nfft:-1:Nfft-GenPar.reb_num+1,:) = H_left;
    end    
        
    Rt = A_sig'*H_ls_comb;
    %Rt = A_sig_us'*H_ls_us;
    Xt = zeros(size(Rt));
    pos_set = [];
        
    % basis search
    for iter_idx = 1:num_iter
        
        if 0
            pdp_est = mean(abs(Rt).^2, 2);
            if 0
            figure(11);
            plot(pdp_est);
            grid on;
            end
            [max_pdp, pos_max] = max(pdp_est);

            if (max_pdp > 1.0*sigma_noise)
                filt_val = (max_pdp - 1.0*sigma_noise) ./ max_pdp;
                %filt_val = 1.0;
            else
                if isempty(pos_set)              
                    filt_val = 1.0;                
                else
                    break;
                end
            end

            pos_set = [pos_set, pos_max];
            delta_X = filt_val.* S(:, pos_max) * Rt(pos_max, :);

            Rt_new = Rt - delta_X;
            %Xt(pos_max,:) = Xt(pos_max,:) + delta_X(pos_max,:);
            Xt(:,:) = Xt(:,:) + delta_X;

            %figure(1); plot(real(Rt(:,1))); hold on; plot(real(Rt_new(:,1))); hold off;
            Rt = Rt_new;
            
            
        end
        
        
        Rtt = S'*Rt;
        pdp_est = mean(abs(Rtt).^2, 2);
        
        [max_pdp, max_pos] = max(pdp_est);
        if max_pdp > sigma_noise
            pos_set = [pos_set, max_pos]
        else
            if isempty(pos_set)
                pos_set = max_pos;
            else
                break;
            end
        end
        
        S_cut = S(:, pos_set);
        X_cut = inv(S_cut'*S_cut)*S_cut' * Rt;
        
        Rt_new = Rt - S_cut * X_cut;
        Rt = Rt_new;
        
        if mean(abs(Rt(:)).^2) < 1.0*sigma_noise
            break;
        end
        norm(Rt_new)
    end
    
    % least square
    pos_unique = unique(pos_set);
    
    A_cut = A_sig(:, pos_unique);   
    A_cut_int = A_sig_int(:, pos_unique);
    
    Xt = inv(A_cut'*A_cut)*A_cut' * H_ls_comb;
    
    pdp_est = mean( abs(Xt).^2 , 2);
    %pdp_est = pdp_est ./ pdp_est;
    W_filt = A_cut_int * inv(A_cut'*A_cut + 0.0*diag(sigma_noise./pdp_est))*A_cut';
    
    H_imp = W_filt * H_ls_comb;
    
end


