% transmit and receive harmonic function on sdr
if exist('plutoTx')
    release(plutoTx)
end

if exist('plutoRx')
    release(plutoRx)
end

close all; clear all; clc;
rng(1);
addpath('../ce');

% sim param
do_tx_only = false;
do_rx_only = true;
do_save_chan = false;
bb_name = 'Dump11.mat'; % Dump1 - bedroom Dump2 far bedroom Dump3 - guestroom sofa, dump 4 - kitchen, dump5-8 balcony (7 on the top)
% dump 9-11 corr
do_load_bb = true;
do_load_chan = 'Dump11.mat'; % 'NearPoint.mat' , 'Dump1.mat'

if do_load_bb
    do_tx_only = false;
    do_rx_only = true;
    do_save_chan = false;
    load(do_load_chan);
end


do_print_ofdm = false;
do_print_ce = true;
num_exp = 5;
num_rep = 1; % number of repeat of data symbols
num_rep_pilot = 1; % number of repeat pilots
do_pdp_print = true;

do_dc_offset = 0;

% tx params
%center_frequency = 3.5e9;
center_frequency = 2.0e9;
sample_rate = 20e6; % 8 * 940e3

db_gain = 0.0;
db_scale = 10^(-db_gain/10);
signal_scale = 1.0 * sqrt(db_scale); % power scale


modu_ord = 4;
num_bits_sym_d = log2(modu_ord);

scen_name_prefix = ['QAM-', num2str(modu_ord), '-scale=', num2str(signal_scale)];

cf_arr = [1];
ce_arr = {'LS', 'IMP',  'HW', 'SW', 'SW-1'};

num_cases = length(cf_arr) * length(ce_arr);

ber_all = cell(1, num_cases);
for i = 1:num_cases
    ber_all{i} = struct;
    ber_all{i}.ber_arr = [];
end

if do_save_chan
    RxBB = cell(num_exp, 1);
end


% OFDM message
delta_t = 1 / sample_rate;

delta_f = 30e3; % subcarrier spacing
frac_guard = 0.5;

N_sc_av = fix( sample_rate / delta_f );
N_fft = 2^( fix(log2(N_sc_av)) );
N_sc_use = fix(N_fft * frac_guard);
guard_length = fix( 0.2 * N_fft );
%guard_length = 0;
CP_len = fix( N_fft * 0.2 );
num_bits_sym = 1;

BlockSize = N_sc_use * num_bits_sym; 
BlockSize_d = N_sc_use * num_bits_sym_d; 

N_repeat = 10; % number of repeats for tx

% create preambule
bark = comm.BarkerCode('SamplesPerFrame',N_fft/2,'Length',13);
preamble = complex(bark());
release(bark);

preamble_full = [preamble; preamble];
preamble_full_cp = [preamble_full(end - CP_len + 1:end); preamble_full];

% pilot frame
data_stream = randi([0, 1], BlockSize, 1);
mod_sym_pilot = 1.0 - 2.0 .* data_stream; % BPSK modulation

tx_ofdm_sym = zeros(N_fft, 1);
tx_ofdm_sym(1+do_dc_offset:N_sc_use/2 + do_dc_offset) = mod_sym_pilot(N_sc_use/2 + 1 : end);
tx_ofdm_sym(end - N_sc_use/2 + 1 : end) = mod_sym_pilot(1:N_sc_use/2);

time_ofdm_sym_pilot = ifft(tx_ofdm_sym) * sqrt(N_fft) * signal_scale;
%time_ofdm_sym_cp_pilot = [time_ofdm_sym_pilot(end - CP_len + 1:end); time_ofdm_sym_pilot];
time_ofdm_sym_pilot_rep = repmat(time_ofdm_sym_pilot, num_rep_pilot, 1);
time_ofdm_sym_cp_pilot = [time_ofdm_sym_pilot(end - CP_len + 1:end); time_ofdm_sym_pilot_rep];

% data frame
data_stream = randi([0, 1], BlockSize_d, 1);
%mod_sym = 1.0 - 2.0 .* data_stream; % BPSK modulation
mod_sym = qammod(data_stream,modu_ord,'InputType','bit','UnitAveragePower',true);

tx_ofdm_sym = zeros(N_fft, 1);
tx_ofdm_sym(1+do_dc_offset:N_sc_use/2 + do_dc_offset) = mod_sym(N_sc_use/2 + 1 : end);
tx_ofdm_sym(end - N_sc_use/2 + 1 : end) = mod_sym(1:N_sc_use/2);

time_ofdm_sym = ifft(tx_ofdm_sym) * sqrt(N_fft) * signal_scale;
%time_ofdm_sym_cp = [time_ofdm_sym(end - CP_len + 1:end); time_ofdm_sym];

time_ofdm_sym_rep = repmat(time_ofdm_sym, num_rep, 1);
time_ofdm_sym_cp = [time_ofdm_sym(end - CP_len + 1:end); time_ofdm_sym_rep];

guard_frame = zeros(guard_length, 1);

full_frame = [preamble_full_cp; time_ofdm_sym_cp_pilot; time_ofdm_sym_cp; guard_frame];
FrameSize = length(full_frame);
sig_len = length(full_frame);

% upsampling and RRC
%sps = sample_rate/data_rate;

%txFilter = comm.RaisedCosineTransmitFilter( ...
%    'OutputSamplesPerSymbol',sps,'FilterSpanInSymbols',2);

%tx_signal = txFilter(IQmessage);

if ~do_rx_only
    figure(11);
    plot(abs(full_frame));
    grid on;

    plutoTx = plutoTxInit(center_frequency,sample_rate);
    transmitRepeat(plutoTx,full_frame);
    % actual Tx

    if (do_tx_only)
        fprintf('Transmission in progress\n');
        keyboard;
        fprintf('Transmission in ended');
        return;
    end
end


ber_cfo_ls = [];
ber_no_cfo_ls = [];

for exp_idx = 1:num_exp
    
    samples_per_frame = 5*sig_len;

    %RTLSDRrx = RTLSDRRxInit(center_frequency,sample_rate,samples_per_frame);
    %rx_signal = RTLSDRrx();
    %release(RTLSDRrx);

    if ~do_load_bb
        plutoRx = plutoRxInit(center_frequency,sample_rate,samples_per_frame);
        rx_signal = plutoRx();
    else
        rx_signal = RxBB{exp_idx};
    end

    if do_save_chan
        RxBB{exp_idx} = rx_signal;
        save(bb_name, 'RxBB');
    end
    %save('rx_sig_m1.mat', 'rx_signal');
    

    if do_print_ofdm
        figure(20);
        plot(10*log10(abs(rx_signal)));
    end

    %frame_rx = frame_rx_cfo(init_pos + 1: init_pos + reception_win);
    frame_rx = rx_signal;
    % step 1 - search begging of the frame
    len_rec = length(frame_rx);
    corr_len = N_fft / 2;

    metric_p = [];

    for start_pos = 0:len_rec-N_fft-1

        A1 = frame_rx(start_pos + 1: start_pos + N_fft/2);
        A2 = frame_rx(start_pos + N_fft/2 + 1: start_pos + N_fft);

        corr = A1'*A2;
        metric_p = [metric_p; corr];

    end

    filt_1 = ones(CP_len, 1);
    metric_p1 = conv(metric_p, filt_1, 'same');


    if do_print_ofdm
        figure(1);
        plot(abs(metric_p));
        grid on

        figure(2);
        plot(abs(metric_p1));
        grid on
    end

    [max_val, pos_start] = max(abs(metric_p1(1:FrameSize)));
    pos_start = pos_start + fix(CP_len / 2);
    %pos_start = CP_len / 2;

    cfo_est = angle( metric_p(pos_start) ) / (2*pi) / (N_fft/2);
    %delta_cfo = cfo_idl - cfo_est;

    x_arr = (0:length(frame_rx)-1).';
    %cfo_est = 0.0;

    k_scen = 1;
    for cf_flag = cf_arr
        
        scen_name = [scen_name_prefix, '-cfo=' , num2str(cf_flag)];
        ber_now = struct;
        
        if (cf_flag == 0)
            cfo_est_use = 0;
        else
            cfo_est_use = cfo_est;
        end
        cfo_corr = exp(1j * (-2*pi*cfo_est_use) * x_arr);
        frame_rx_cfo = frame_rx .* cfo_corr;


        frame_ta = frame_rx_cfo( pos_start : pos_start + N_fft + (CP_len + num_rep_pilot * N_fft) + (CP_len + num_rep*N_fft) );

        % compensate CFO
        %x_arr = (0:length(frame_ta)-1).';
        %cfo_corr = exp(1j * (-cfo_est) * x_arr);
        %frame_ta = frame_ta .* cfo_corr;

        start_pos = N_fft + CP_len + 1;
        end_pos = start_pos - 1 + num_rep_pilot * N_fft;
        pilot_frame = frame_ta(start_pos : end_pos);
        
        % combining of pilots
        pilot_frame = reshape(pilot_frame, N_fft, num_rep_pilot);
        pilot_frame = mean(pilot_frame, 2);
        
        pilot_freq = fft(pilot_frame);


        pilot_freq_bb = [pilot_freq(end - N_sc_use/2 + 1 : end) ; pilot_freq(1+do_dc_offset : N_sc_use/2 + do_dc_offset)];
        
        % noise outside the band
        pilot_noise = pilot_freq(N_sc_use/2 + 2 : end - N_sc_use/2);
        sigma_out_band = norm(pilot_noise(:))^2 / length(pilot_noise(:));
        
        % rough SNR est
        sigma_es_band = norm(pilot_freq_bb(:))^2 / length(pilot_freq_bb(:));
        
        if sigma_es_band > sigma_out_band
            SNR_sc_rough = 10*log10( (sigma_es_band - sigma_out_band)/sigma_out_band);
            fprintf('expIdx=%d SNR=%f\n', exp_idx, SNR_sc_rough);
        else
            fprintf('expIdx=%d Cant estimate SNR Es < N0\n', exp_idx);
        end
        
        %mod_sym_pilot

        % Least square CE
        % Least square CE
        scen_name_partial = scen_name;
        for ce_method = ce_arr
            
            h_ls = pilot_freq_bb ./ mod_sym_pilot;

            if do_print_ce && strcmp(ce_method{1}, 'LS')
                figure(420 + exp_idx);
                subplot(2,1,1);
                plot(real(h_ls));
                legend('Real LS');
                grid on;
                subplot(2,1,2);
                plot(imag(h_ls));
                legend('Imag LS');
                grid on;
            end


            % HW/SW CE
            scen_name = [scen_name_partial, '-CE=', ce_method{1}];
            
            % INIT CE params
            N_sc_ls = length(h_ls);
            N_fft_ce = 2^( ceil(log2(N_sc_ls) ) + 1 );
            
            
            % system parameters
            GenPar.Nsc = length(h_ls);
            GenPar.comb = 1;
            GenPar.delta_f = delta_f; % supcarrier spacing
            GenPar.Nrx = 1;
            GenPar.SNR = -8;
            GenPar.do_quad = 0;
            GenPar.do_pow2 = true;
            %GenPar.Nexp_total = 30;

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
            ChanInfo.Nfft = N_fft_ce;
            %ChanInfo.win_max = fix( 1*ChanInfo.tau_max * GenPar.delta_f * ChanInfo.Nfft  );
            ChanInfo.win_max = 70;
            %ChanInfo.win_guard = 12; % max sync error for Nfft=2048
            %ChanInfo.win_min = fix( ChanInfo.win_guard * ChanInfo.Nfft / 2048);
            ChanInfo.win_min = 50;
            % end INIT
            
            if ~strcmp(ce_method{1}, 'LS')
            
                % Channel estimation
                win_for = 70;
                win_back = 50;

                
                
                h_ls_zp = zeros(N_fft_ce, 1);
                h_ls_zp(1:N_sc_ls) = h_ls;

                % rebuilding
                reb_right = 6;
                reb_left = 6;


                ce_name = ce_method{1};


                if strcmp(ce_method{1}, 'HW-Reb') || strcmp(ce_method{1}, 'SW-Reb') 
                    do_reb = true;
                else
                    do_reb = false;
                end


                if do_reb
                    % right
                    h_ls_zp(N_sc_ls + 1 : N_sc_ls + reb_right) = h_ls_zp(N_sc_ls);

                    %left
                    h_ls_zp(end - reb_left + 1 : end) = h_ls_zp(1);
                end

                h_time = ifft(h_ls_zp) * sqrt(N_fft_ce);

                % noise est;
                h_noise = h_time(win_for+1:end-win_back);
                sigma_noise = mean(abs(h_noise).^2);

                filt_hw = zeros(N_fft_ce, 1);
                filt_hw(1:win_for) = 1.0;
                filt_hw(N_fft_ce - win_back + 1 : end) = 1.0;

                h_time_filt = h_time .* filt_hw;

                % soft win
                pdp_time = abs(h_time_filt).^2;
                pos_high = find(pdp_time  > sigma_noise);
                filt_sw = zeros(N_fft_ce, 1);
                filt_sw(pos_high) = ( pdp_time(pos_high) - sigma_noise) ./ pdp_time(pos_high);

                % sw
                if strcmp(ce_method{1}, 'SW') || strcmp(ce_method{1}, 'SW-Reb')
                    h_time_filt = h_time_filt .* filt_sw;
                end

                h_freq_us = fft(h_time_filt) ./ sqrt(N_fft_ce);   
                
                
                % IMP method
                if strcmp(ce_method{1}, 'IMP')
                    [H_imp, Info] = imp_est(h_ls, GenPar, ChanInfo);
                    h_ls = H_imp(1:N_sc_ls);
                elseif strcmp(ce_method{1}, 'SW-1')
                    [H_sw, Info] = sw_est(h_ls, GenPar, ChanInfo);
                    h_ls = H_sw(1:N_sc_ls);
                else
                    h_ls = h_freq_us(1:N_sc_ls);
                end
                
                
                
                % chane est end
            end
            
            h_time = ifft(h_ls, N_fft_ce);
            pdp_now = abs(h_time).^2;
            
            [pdp_sort] = sort(pdp_now);
            len_corr = 8;
            h_mat = reshape(h_ls, len_corr, []);
            Rhh = h_mat*h_mat';
            [V,S,~] = svd(Rhh);
            r1 = S(1,1) / S(2,2);
                
            if do_pdp_print && strcmp(ce_method{1}, 'SW-1') && cf_flag == 1
                figure(k_scen + exp_idx * 10);
                w_1 = 10;
                w_2 = 30;
                plot([pdp_now(end-w_1+1:end); pdp_now(1:w_2)]);
                grid on;
                legend(ce_method{1});
                %pause(1);
                %fprintf('Pause in reception');
            end
            
            
            % Post CE SNRest
            %h_ls = pilot_freq_bb ./ mod_sym_pilot;
            if ~strcmp(ce_method{1}, 'LS')
                u_samples = pilot_freq_bb - h_ls .* mod_sym_pilot;
                
                SNR_ce = 10*log10( norm(h_ls(:))^2 / norm(u_samples(:))^2 );
                fprintf('expIdx=%d postCE=%s SNR=%f\n', exp_idx, ce_method{1}, SNR_ce);
            end
            

            % data frame
            start_pos = (N_fft + CP_len) + (num_rep_pilot * N_fft + CP_len + 1);
            end_pos = start_pos + num_rep * N_fft - 1;
            data_frame = frame_ta(start_pos : end_pos);
            
            % combining
            data_frame = reshape(data_frame, N_fft, num_rep);
            data_frame = mean(data_frame, 2);
            
            data_freq = fft(data_frame);


            data_freq_bb = [data_freq(end - N_sc_use/2 + 1 : end) ; data_freq(1+do_dc_offset : N_sc_use/2 + do_dc_offset)];

            %mod_sym

            demod_sym = data_freq_bb ./ h_ls;

            %demod_sym = [demod_sym(N_sc_use/2 + 1 : end); demod_sym(1:N_sc_use/2)];

            evm_est = norm(demod_sym - mod_sym) ./ norm(mod_sym);

            if strcmp(ce_method{1}, 'HW') && 0
                figure(6 + cf_flag*10);
                scatterplot(demod_sym);
                grid on;
                title(['CFOcomp=', num2str(cf_flag)])
            end

            %bit_dec = double(demod_sym < 0); % BPSK demodulation
            bit_dec = qamdemod(demod_sym,modu_ord, 'UnitAveragePower', true,  'OutputType','bit');
            ber = sum(bit_dec ~= data_stream) / length(data_stream);

            
            
            ber_all{k_scen}.ber_arr = [ber_all{k_scen}.ber_arr, ber];
            ber_all{k_scen}.scen_name = scen_name;
            
            %fprintf('Current : ExpIdx=%d Scen=%s BER=%f EVMest=%f\n', exp_idx, scen_name, ber, evm_est);
            
            k_scen = k_scen + 1;
        end
    end

    for scen_idx = 1:num_cases        
        fprintf('Num avg=%d Name=%s BER=%f\n', exp_idx, ber_all{scen_idx}.scen_name, mean( ber_all{scen_idx}.ber_arr ));
    end
end

if exist('plutoRx')
    release(plutoRx);
end
