function [H_left, H_right] = calc_rebuild(H_ls_comb, GenPar, ChanInfo)
% calc rebuilding

    [Nsc_comb, Nrx] = size(H_ls_comb);

    reb_num = GenPar.reb_num;
    reb_base = GenPar.reb_base;
    
    x_arr_reb = [-1:-1:-reb_num];
    x_arr_base = [0:1:reb_base-1];
    
    % left
    [R_self] = calc_cov_wien_mx(x_arr_base, x_arr_base, GenPar, ChanInfo);
    [R_cross] = calc_cov_wien_mx(x_arr_reb, x_arr_base, GenPar, ChanInfo);
    
    linInvSNR = 10^(-GenPar.SNR/10);
    W_filt =  R_cross' * pinv(R_self + eye(reb_base)*linInvSNR);
    
    H_base = H_ls_comb(1:reb_base, :);
    H_left = W_filt * H_base;
    
    % right
    x_arr_base = Nsc_comb-reb_base+1:Nsc_comb;
    x_arr_reb = Nsc_comb+1:Nsc_comb+reb_num;
    
    [R_self] = calc_cov_wien_mx(x_arr_base, x_arr_base, GenPar, ChanInfo);
    [R_cross] = calc_cov_wien_mx(x_arr_reb, x_arr_base, GenPar, ChanInfo);
    
    linInvSNR = 10^(-GenPar.SNR/10);
    W_filt = R_cross' * pinv(R_self + eye(reb_base)*linInvSNR);
    
    H_base = H_ls_comb(x_arr_base, :);
    H_right = W_filt * H_base;
    
    
end

function [R_out] = calc_cov_wien_mx(x_arr_reb, x_arr_base, GenPar, ChanInfo)
    
    T1 = repmat(x_arr_base(:), 1, length(x_arr_reb));
    T2 = repmat( x_arr_reb, length(x_arr_base), 1);
    D = T2 - T1;
    
    R_out = 1.0 ./ (1 - 1.0*2*pi*1j*(D*GenPar.delta_f*GenPar.comb*(ChanInfo.tau_max/3.0)));

end