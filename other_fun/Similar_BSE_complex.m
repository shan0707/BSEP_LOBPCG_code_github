function [X1, X2, Lambda] = Similar_BSE_complex(A, B, k_pos)
%% compute H = [A, B; -conj(B), -conj(A)], where A = A', B = -B.'
%% X = [X1, conj(X2); 
%%      X2, conj(X1)]
%% HX = Xdiag(Lambda, -Lambda), where X'X=I
% output k_pos largest
i = sqrt(-1);
n = size(A, 1);

At = [imag(A), -real(A); 
      real(A),  imag(A)];
Bt = [imag(B),  real(B); 
      real(B), -imag(B)];
M = At + Bt;
M = (M - M')/2;

[U, T] = hess(M);
t = diag(T, 1);
T = diag(t, 1) + diag(t, -1); % = sqrt(-1)*D1'*T*D1 such that T is symmetric
[V, D] = eig(T);

if nargin > 2
    index = 2*n : -1 : 2*n - (k_pos - 1);
else
    index = 2*n : -1 : n + 1;
end
Lambda = D(index, index);
V = V(:, index);

d1 = (-i).^(1 : 2*n).';

%X = Q*U*D*V
S = U*(d1.*V);
S1 = S(1 : n, :);
S2 = S(n + 1 : 2*n, :);
S2i = i*S2;
X1 = (S1 - S2i)/sqrt(2);
X2 = (S1 + S2i)/sqrt(2);