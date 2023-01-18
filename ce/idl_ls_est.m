function [H_idl_ls, Info] = idl_ls_est(H_noisy, GenPar, ChanInfo)
% ideal LS denoise
    H_ls_comb = H_noisy(1:GenPar.comb:end,:);
    
    comb = GenPar.comb;
    Nsc_comb = GenPar.Nsc / GenPar.comb;
    
    freq_idl = ChanInfo.freq_idl';
    freq_idl = freq_idl .* comb;
    
    x_arr = 0:Nsc_comb-1;
    x_arr = x_arr(:);
    
    x_int = 0:(comb*Nsc_comb-1);
    x_int = x_int ./ comb;
    x_int = x_int(:);
    
    sigma_noise_sc = GenPar.sigma_noise_sc;
    
    A_idl = exp(-1j*2*pi*x_arr*freq_idl) ./ sqrt(Nsc_comb);
    A_idl_int = exp(-1j*2*pi*x_int * freq_idl) ./ sqrt(Nsc_comb);
    
    h_time = pinv(A_idl'*A_idl)*A_idl'*H_ls_comb;
    pdp = mean(abs(h_time).^2 , 2);
    %pdp = pdp ./ abs(pdp);
    invPDP = sigma_noise_sc ./ pdp;
    
    %H_idl_ls = A_idl_int * pinv(A_idl'*A_idl + 1.0*diag(invPDP))*A_idl'*H_ls_comb;   
    
    % blockwise
    Nbl = 48*2/GenPar.comb;
    num_gr = Nsc_comb / Nbl;
    A_idl_cut = A_idl(1:Nbl, :);
    A_idl_cut_int = A_idl_int(1:Nbl*comb, :);
    W_cut = A_idl_cut_int * pinv(A_idl_cut'*A_idl_cut + 1.0*diag(invPDP))*A_idl_cut';
    for gr_idx = 1:num_gr
        H_cut = H_ls_comb( (gr_idx - 1)*Nbl + 1: gr_idx * Nbl, : );
        
        H_idl_ls((gr_idx - 1)*Nbl*comb + 1: gr_idx * Nbl*comb, : ) = W_cut * H_cut;        
    end
   
    
    
    
    
    H_tmp = A_idl*pinv(A_idl'*A_idl + 1.0*diag(invPDP))*A_idl'*H_ls_comb;
    
    A_dft = dftmtx(Nsc_comb) ./ sqrt(Nsc_comb);
    A_dft_int = dftmtx(comb*Nsc_comb) ./ sqrt(Nsc_comb);
    
    h1 = A_dft' * H_tmp;
    h2 = zeros(size(h1,1)*comb, size(h1,2));
    
    h2(1:Nsc_comb/2, :) = h1(1:Nsc_comb/2, :);
    h2(end - Nsc_comb/2 + 1:end, :) = h1(end - Nsc_comb/2 + 1:end, :);
    
    h3 = A_dft_int * h2;
    
    %H_idl_ls = h3;
    
    Info = [];
end

