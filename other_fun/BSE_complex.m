function [X1, X2, Lambda] = BSE_complex(A, B, k_pos)
%% compute BSE eigenproblem
i = sqrt(-1);
n = size(A, 1);

At = [real(A), imag(A); 
     -imag(A), real(A)];
Bt = [real(B), -imag(B); 
     -imag(B), -real(B)];
M = At + Bt;
M = (M + M')/2;

[L, flag] = chol(M, 'lower');
if flag ~= 0
    if any(~isfinite(M), 'all')
        error('BSE_complex:NonFiniteProjection', ...
            'The projected BSE matrix contains NaN or Inf values.');
    end
    lambda_min = min(real(eig(M)));
    scale = max(1, norm(M, 1));
    shift = max(1e-10, -lambda_min + 1e-12*scale);
    identity = eye(2*n, 'like', M);
    [L, flag] = chol(M + shift*identity, 'lower');
    attempts = 0;
    while flag ~= 0 && attempts < 8
        shift = 10*shift;
        [L, flag] = chol(M + shift*identity, 'lower');
        attempts = attempts + 1;
    end
    if flag ~= 0
        error('BSE_complex:CholeskyFailure', ...
            'Projected matrix is not positive definite after shift %.3e.', shift);
    end
end

W = L'*Jmul(L);
W = (W - W')/2;
[U, T] = hess(W);
t = -diag(T, -1);
T = diag(t, 1) + diag(t, -1);
[V, D_eig] = eig(T);

if nargin > 2
    index = n + 1 : n + k_pos;
else
    index = n + 1 : 2*n;
end

Lambda = D_eig(index, index);

d1 = i.^(0 : 2*n - 1).'; 
%X = Q*(L'\U)*(d1.*V(:, index));
S = (L'\U)*(d1.*V(:, index));
S1 = S(1 : n, :);
S2 = S(n + 1 : 2*n, :);
sqrt_lambda = diag(Lambda).^(1/2);
X1 = (S1 - i*S2)/sqrt(2).*sqrt_lambda.';
X2 = (S1 + i*S2)/sqrt(2).*sqrt_lambda.';
end
