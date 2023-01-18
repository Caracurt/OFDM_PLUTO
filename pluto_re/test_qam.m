% test QAM mod/demod

modu_ord = 16;
num_bits = log2(modu_ord);

data = randi([0 1],1000*num_bits,1);
txSig = qammod(data,modu_ord,'InputType','bit','UnitAveragePower',true);

scatterplot(txSig)
title('64-QAM, Average Power = 1 W')

z = qamdemod(txSig,modu_ord, 'UnitAveragePower', true,  'OutputType','bit');
s = isequal(data, double(z));