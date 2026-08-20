function [X, Y] = mysqrtm(A)
%% A, X, and Y are Hermitian positive definite and of the form [a, b; conj(b), a] 
% such that A = X^2 = Y^(-2),  x are real
a = A(1, 1);
b = A(1, 2);
r = abs(b);
u = sqrt(a + r);
v = sqrt(a - r);
x = (u + v)/2;
y = b/(u + v);
X = [x, y; conj(y), x];
Y = [x, -y; -conj(y), x]/(u*v);