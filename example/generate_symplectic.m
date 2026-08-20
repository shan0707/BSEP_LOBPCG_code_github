function M = generate_symplectic(n)
% parameters for symplectic Gauss transformation
Gauss.c = 1.2;
Gauss.m = round(n/5);
Gauss.d = -sqrt(Gauss.m);
DD = diag([1 : n 1 : n]);
rng default
Q1 = randn(n, n);
Q2 = randn(n, n);
[U, ~, ~] = svd(Q1 + sqrt(-1)*Q2);
M = [real(U) -imag(U); imag(U) real(U)];
d1 = [ones(1, Gauss.m - 2) Gauss.c Gauss.c ones(1, n - Gauss.m)];
d2 = [ones(1, Gauss.m - 2) 1/Gauss.c 1/Gauss.c ones(1, n - Gauss.m)];
L2 = zeros(n, n); L2(Gauss.m, Gauss.m-1) = Gauss.d;
L2(Gauss.m - 1, Gauss.m) = Gauss.d;
Lmcd = [diag(d1) L2; zeros(n, n) diag(d2)];%sympl. Gauss transformation type I
M = M*Lmcd;
M = 0.5*((M*DD)*M' + M*(DD*M'));
M = 0.5*(M + M');