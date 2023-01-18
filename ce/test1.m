% test
Nrx = 64;
Nl = 12;
H = randn(Nrx, Nl) + 1j*randn(Nrx, Nl);

exp_num = 100*12*128;
for exp_idx = 1:exp_num
    H = randn(Nrx, Nl) + 1j*randn(Nrx, Nl);
    [Q,R] = qr(H);
end