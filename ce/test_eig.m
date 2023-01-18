% test eig dist of noise

Nrx = 64;
Nsc = 96*1;
Nexp = 100;
SNR = 30;
sigma0 = 10^(-SNR/10);

eta1_arr = [];
eta2_arr = [];
eta3_arr = [];

Nbase = 12;
rho_arr = [];

IOT_sigma = 10;

rng(1);
for exp_idx = 1:Nexp
    
    H_idl = randn(Nrx, 1) + 1j*randn(Nrx, 1);
    
    V = randn(Nrx, Nbase) + 1j*randn(Nrx, Nbase);
    [Q1,R1] = qr(V, 0);
    
    H = sqrt(sigma0/2) .* ( randn(Nrx, Nsc) + 1j*randn(Nrx, Nsc) );

    alpha = (IOT_sigma./sqrt(2)) .* ( randn(Nbase, Nsc) + 1j*randn(Nbase, Nsc)  );
    
    U = Q1 * alpha;
    
    U_full = U + H;
    
    R = U_full*U_full' ./ Nsc;

    [U,S,V] = svd(R);

    l = diag(S);

    eta1 = sigma0 / l(end);
    eta1_arr = [eta1_arr, eta1];

    eta2 = l(1) / l(end);
    
    eta3 = l(end) * 10 / sigma0;
    eta3_arr = [eta3_arr, eta3];
    
    eta2_arr = [eta2_arr, eta2];
           
    % receiver
    load = max(abs(diag(R)))*2^(-12);
    
    R_diag = R + eye(Nrx)*load;
    w = inv(R_diag)*H_idl;
    
    rho = abs(w'*H_idl)^2;
    
    
    rho_arr = [rho_arr, rho];
    
end

mean(rho_arr)

figure(1);
plot(eta1_arr);
grid on;

figure(2);
plot(eta2_arr);
grid on;

figure(3);
plot(eta3_arr);
grid on;



