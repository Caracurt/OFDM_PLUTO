% test quadratic regression

a_vect = 3*randn(3, 1);

a = -abs(a_vect(1)); b = 3*a_vect(2); c = a_vect(3) + 3;

Sm1 = a - b + c;
S = c;
Sp1 = a + b + c;

b_est = (Sp1 - Sm1)/2;
a_est = (Sp1 + Sm1- 2*S)/2;

delta = (Sp1 - Sm1) / ( 2.0*(2*S - Sp1 - Sm1) );
delta 
delta_idl = -b/(2*a)
x_arr = -10:0.1:10;

y_arr = a.*x_arr.^2 + b * x_arr + c;
max(y_arr)
figure(1);
plot(x_arr, y_arr, '-r');
grid on;

y_delta = a * delta^2 + b * delta + c
y_delta_idl = a * delta_idl^2 + b * delta_idl + c