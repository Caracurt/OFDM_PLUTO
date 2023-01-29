function [H_est, Info] = est_est(h_ls, EspPar)
% ESPRIT channel estimation

    d = EspPar.d;
    Np = EspPar.Np;

    Nsc = length(h_ls);


    % collect cov matrix
    Y = [];
    for s_idx = 1:Nsc - Np
        Y_c = h_ls( s_idx : s_idx + Np - 1 );
        Y = [Y, Y_c];
    end
    Ryy = Y*Y';

    [U,S,V] = svd(Ryy);

    U_up = U(1:end-1,1:d);
    U_low = U(2:end,1:d);

    T = pinv(U_up'*U_up)*U_up'*U_low;
    [Ui, lambd] = eig(T);
    lambd = diag(lambd);

    tau = angle(lambd);
    tau = reshape(tau, 1, d);

    x_arr = 0:Nsc-1;
    x_arr = x_arr.';

    Exp_arr = exp(1j*x_arr*tau);

    W = pinv(Exp_arr'*Exp_arr) * Exp_arr' * h_ls;

    H_est = Exp_arr * W;
    Info = [];

end