% testing platform for channel estimation
clear all;


% system parameters
GenPar.Nsc = 32*12;
GenPar.comb = 2;
GenPar.delta_f = 30e3; % supcarrier spacing
GenPar.Nrx = 16;
GenPar.SNR = -8;
GenPar.do_quad = 0;
GenPar.do_pow2 = true;
GenPar.Nexp_total = 30;

% channel parameters
ChanInfo.num_taps = 6;
ChanInfo.tau_max = 5.5e-6;

GenPar.do_rebuild = true;
GenPar.reb_num = 6;
GenPar.reb_base = 12;

% OMP params
GenPar.us_factor_omp = 1.0;
GenPar.num_omp_iter = 12;

% ISTA params
GenPar.ista_iter = 12;
GenPar.do_amp = 0;
GenPar.ista_iter_warm = 1;


Nfft = GenPar.us_factor_omp * max(64, 2^(ceil(log2(GenPar.Nsc/GenPar.comb)))  );
ChanInfo.Nfft = Nfft;
ChanInfo.win_max = fix( 1*ChanInfo.tau_max * GenPar.delta_f * ChanInfo.Nfft  );
ChanInfo.win_guard = 12; % max sync error for Nfft=2048
ChanInfo.win_min = fix( ChanInfo.win_guard * ChanInfo.Nfft / 2048);



%% CHANEL GENRATION
rng(3);
nmse_arr = [];
nmse_idl = [];
nmse_arr_omp = [];
nmse_arr_ista = [];
nmse_arr_imp = [];
nmse_arr_timp = [];

for exp_idx = 1:GenPar.Nexp_total
    GenPar.exp_idx = exp_idx;
    [H_idl, Info] = gen_channel(GenPar, ChanInfo);

    %% ADD NOISE
    sigma_noise = 10^(-GenPar.SNR/10);
    noise_vect = sqrt(sigma_noise/2) .* (randn(GenPar.Nsc, GenPar.Nrx) + 1j*randn(GenPar.Nsc, GenPar.Nrx));
    noise_vect = sqrt(sigma_noise) * noise_vect ./ norm(noise_vect);

    sigma_noise_sc = norm(noise_vect(:))^2 ./ length(noise_vect(:));
    GenPar.sigma_noise_sc = sigma_noise_sc;
    
    H_noisy = H_idl + noise_vect;
    
    % SNR est
    SNR_est = 10*log10( norm(H_idl(:))^2 / (norm(H_idl(:))^2 + norm(noise_vect(:))^2) );
    

    %% MAKE ESTIMATION

    % idela freq denoise (bound
    ChanInfo.freq_idl = Info.freq_out;
    
    [H_idl_ls, Info] = idl_ls_est(H_noisy, GenPar, ChanInfo);
    nmse_idl_c = sqrt( norm(H_idl_ls(:) - H_idl(:))^2 / norm(H_idl(:))^2 );
    nmse_idl = [nmse_idl, nmse_idl_c];
    
    % soft window estimation
    [H_sw, Info] = sw_est(H_noisy, GenPar, ChanInfo);

    %[H_sw, Info] = dctw_est(H_noisy, GenPar, ChanInfo);
    
    nmse_sw_c = sqrt( norm(H_sw(:) - H_idl(:))^2 / norm(H_idl(:))^2 );
    
    nmse_arr = [nmse_arr, nmse_sw_c];
    
    % OMP channel estimate
    [H_omp, Info] = omp_est(H_noisy, GenPar, ChanInfo);
    
    nmse_omp_c = sqrt( norm(H_omp(:) - H_idl(:))^2 / norm(H_idl(:))^2 );
    
    nmse_arr_omp = [nmse_arr_omp, nmse_omp_c];
    
    % ISTA
    [H_ista, Info] = ista_est(H_noisy, GenPar, ChanInfo);
    
    nmse_ista_c = sqrt( norm(H_ista(:) - H_idl(:))^2 / norm(H_idl(:))^2 );    
    nmse_arr_ista = [nmse_arr_ista, nmse_ista_c];
    
    % IMP 
    [H_imp, Info] = imp_est(H_noisy, GenPar, ChanInfo);
    
    H_imp_use = H_imp(12:end-12,:);
    H_idl_use = H_idl(12:end-12,:);
    nmse_imp_c = sqrt( norm(H_imp_use(:) - H_idl_use(:))^2 / norm(H_idl_use(:))^2 );    
    nmse_arr_imp = [nmse_arr_imp, nmse_imp_c];
    
    % TIMP
    [H_timp, Info] = timp_est(H_noisy, GenPar, ChanInfo);
    
    nmse_imp_c = sqrt( norm(H_timp(:) - H_idl(:))^2 / norm(H_idl(:))^2 );    
    nmse_arr_timp = [nmse_arr_timp, nmse_imp_c];
    
    
    % debug
    if exp_idx == 1
        figure(100);
        plot(abs(ifft(H_noisy(:,1))));
        ff=1;
        figure(1);
        x_arr = 0:GenPar.Nsc-1;
        plot(x_arr(:), real(H_idl(:,1)), 'r', 'LineWidth', 2.0);
        hold on;
        plot(x_arr(:), real(H_noisy(:,1)), 'k', 'LineWidth', 2.0);
        hold on
        plot(x_arr(:), real(H_sw(:,1)), 'b', 'LineWidth', 2.0);
        grid on;
        plot(x_arr(:), real(H_omp(:,1)), 'm', 'LineWidth', 2.0);
        grid on;
        plot(x_arr(:), real(H_ista(:,1)), 'g', 'LineWidth', 2.0);
        grid on;
        plot(x_arr(:), real(H_imp(:,1)), 'y', 'LineWidth', 2.0);
        grid on;
        plot(x_arr(:), real(H_timp(:,1)), '-sk', 'LineWidth', 2.0);
        grid on;
        plot(x_arr(:), real(H_idl_ls(:,1)), '-sb', 'LineWidth', 2.0);
        grid on;
        legend('IDL' , 'LS', 'SW', 'OMP', 'ISTA', 'IMP', 'TIMP', 'IDLLS');
        hold off;
    end
end

fprintf('Average NMSE SW=%f\n', mean(nmse_arr));
fprintf('Average NMSE OMP=%f\n', mean(nmse_arr_omp));
fprintf('Average NMSE ISTA=%f\n', mean(nmse_arr_ista));
fprintf('Average NMSE IMP=%f\n', mean(nmse_arr_imp));
fprintf('Average NMSE TIMP=%f\n', mean(nmse_arr_timp));
fprintf('Average NMSE IDLLS=%f\n', mean(nmse_idl));
