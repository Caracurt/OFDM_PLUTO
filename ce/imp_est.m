function [H_imp, Info] = imp_est(H_noisy, GenPar, ChanInfo)
% imp channel estimation
% for SDR test iter=14
    Info = [];    
    H_ls_comb = H_noisy(1:GenPar.comb:end,:);
    
    do_dct = 0;
    
    comb = GenPar.comb;
    Nsc_comb = GenPar.Nsc / GenPar.comb;
    Nfft = ChanInfo.Nfft;

    num_iter = GenPar.ista_iter;
    
    A_dft = dftmtx(Nfft) ./ sqrt(Nfft);
    A_dft_int = dftmtx(Nfft*GenPar.comb) ./ sqrt(Nfft);
        
    
    % wild try
    if do_dct
        %A_dct = dct(eye(Nfft));
        A_dct = real(A_dft) * sqrt(2);
        %A_dct_int = dct(Nfft*GenPar.comb) ./ sqrt(comb);
        A_dct_int = real(A_dft_int) * sqrt(2);
    end
    
    % init determine window sizes
    %win_max = fix( ChanInfo.tau_max * GenPar.delta_f * Nfft * GenPar.comb  );
    win_max = ChanInfo.win_max;
    %win_min = fix( ChanInfo.win_guard * Nfft * GenPar.comb / 512);
    win_min = ChanInfo.win_min;
    
    %pos_signal = [(Nfft - win_min + 1:Nfft), (1:win_max)];
    pos_signal = [(1:win_max), (Nfft - win_min + 1:Nfft)];
    
    %freq_sig = pos_signal ./ Nfft;
    pos_sig_int = pos_signal;
    
    pos_high = find(pos_sig_int > Nfft/2);
    pos_sig_int(pos_high) = pos_sig_int(pos_high) - Nfft;
    
    freq_sig = (pos_sig_int-1) ./ Nfft;
    %pos_signal = pos_sig_int;
    
    x_int = (0:GenPar.Nsc-1).';
    x_int = x_int ./ comb;
    A_sig_int = exp( -1j * 2*pi * x_int * (pos_sig_int-1) ./ Nfft ) ./ sqrt(Nfft);
    
    % temporal
    %A_sig_int = A_dft(1:GenPar.Nsc, pos_signal);
    
    
    A_sig = A_dft(1:Nsc_comb, pos_signal);
    
    if do_dct
        A_sig_use = A_dct(1:Nsc_comb, pos_signal);
    else
        A_sig_use = A_sig;
    end
    
    A_sig_us = A_dft(:, pos_signal);
    
    win_all = length(pos_signal);
    
    pos_noise = [win_max+1:Nfft-win_min];
    A_noise = A_dft(1:Nsc_comb, pos_noise);
    
    sigma_est_noise = A_noise'*H_ls_comb;
    sigma_est_noise = mean( mean(abs(sigma_est_noise).^2, 2) );
    
    sigma_noise = sigma_est_noise;
    
    S = A_sig_use'*A_sig_use;
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
        
    Rt = A_sig_use'*H_ls_comb;
    %Rt = A_sig_us'*H_ls_us;
    Xt = zeros(size(Rt));
    pos_set = [];
    pos_set_quad = [];
        
    % basis search
    for iter_idx = 1:num_iter
        
        pdp_est = mean(abs(Rt).^2, 2);
        if 0
        figure(11);
        plot(pdp_est);
        grid on;
        end
        [max_pdp, pos_max] = max(pdp_est);
        
        
        not_add = 0;
        if (max_pdp > 1.0*sigma_noise)
            
            if sum(find(pos_set == pos_max)) > 0
                not_add = 1;
            end
            
            filt_val = (max_pdp - 1.0*sigma_noise) ./ max_pdp;
            %filt_val = 1.0;
            
            % parabolic regression
            if (pos_max ~= 1) && (pos_max ~= win_all) && GenPar.do_quad
                
                S0 = pdp_est(pos_max);
                Sp1 = pdp_est(pos_max + 1);
                Sm1 = pdp_est(pos_max - 1);
                
                delta = (Sp1 - Sm1) / ( 2.0*(2*S0 - Sp1 - Sm1) );
                
                freq_max = freq_sig(pos_max) + (0.5/Nfft)*delta;
                %freq_max = freq_sig(pos_max) + (1.0/Nfft)*delta;
                
                a_est = (Sp1 + Sm1 - 2*S0) / 2;
                b_est = (Sp1 - Sm1)/2;
                c_est = S0;
                
                
                pdp_max_quad = a_est * delta^2 + b_est*delta + c_est;
                
                %filt_val = (pdp_max_quad - 1.0*sigma_noise)/pdp_max_quad;
                %filt_val = 1.0;
                %pos_set_quad = [pos_set_quad, freq_max];
                
                a_basis = exp(-1j*2*pi*freq_max*(0:Nsc_comb-1).') ./ sqrt(Nsc_comb);                
                s_cancel = A_sig' * a_basis * sqrt(Nfft/Nsc_comb);
                
                %s_cancel = S(:,pos_max);
            else
                s_cancel = S(:,pos_max);
                freq_max = freq_sig(pos_max);                
            end
            
            %filt_val = 1.0;
        else
            if isempty(pos_set)              
                filt_val = 1.0;  
                freq_max = freq_sig(pos_max);
            else
                break;
            end
        end
        
        if not_add == 0
            pos_set_quad = [pos_set_quad, freq_max];
            pos_set = [pos_set, pos_max];
        end
        
        %delta_X = filt_val.* S(:, pos_max) * Rt(pos_max, :);
        delta_X = filt_val.* s_cancel * Rt(pos_max, :);
        
        Rt_new = Rt - delta_X;
        %Xt(pos_max,:) = Xt(pos_max,:) + delta_X(pos_max,:);
        Xt(:,:) = Xt(:,:) + delta_X;
        
        %figure(1); plot(real(Rt(:,1))); hold on; plot(real(Rt_new(:,1))); hold off;
        Rt = Rt_new;
    end
    
    % least square
    %pos_unique = unique(pos_set);
    pos_unique = (pos_set);
    
    A_cut = A_sig(:, pos_unique);   
    A_cut_int = A_sig_int(:, pos_unique);
    
    pdp_est = mean( abs(Xt(pos_unique,:)).^2 , 2);
    %pdp_est = pdp_est ./ pdp_est;
    %W_filt = A_cut_int * inv(A_cut'*A_cut + 0.75*diag(sigma_noise./pdp_est))*A_cut';
    W_filt = A_cut_int * inv(A_cut'*A_cut + 3.0*diag(sigma_noise./pdp_est))*A_cut';
    
    H_imp = W_filt * H_ls_comb;
    
    
    pos_set_quad = unique(pos_set_quad);
    
    A_basis = exp(-1j*2*pi*(0:Nsc_comb-1).' * pos_set_quad) ./ sqrt(Nsc_comb);
    
    x_arr_comb = (0:comb*Nsc_comb-1)./(comb);
    x_arr_comb = x_arr_comb(:);
    A_basis_int = exp(-1j*2*pi*x_arr_comb * pos_set_quad) ./ sqrt(Nsc_comb);
    
    if GenPar.do_quad
        W_filt1 = A_basis_int * inv(A_basis'*A_basis + 0.3*diag(sigma_noise./pdp_est))*A_basis';    
        H_imp = W_filt1 * H_ls_comb;
    end

