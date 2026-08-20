function [UX, UY] = ortho_C_SVQB(UX, UY)
%% Cn-orthonormalize [UX, conj(UY); UY, conj(UX)] using SVQB
k = size(UX, 2);
conj_UX = conj(UX);
conj_UY = conj(UY);
A = UX'*UX - UY'*UY;
T = UX'*conj_UY;
B = T - T.';

diagA = real(diag(A));
scale_floor = eps(max(1, max(abs(diagA))));
d = 1./sqrt(max(abs(diagA), scale_floor));
A = d.*A.*d.';
B = d.*B.*d.';
A = (A + A')/2;
B = (B - B.')/2;

%[VX, VY, Lambda] = Similar_BSE_complex(A, B);
%V = [VX; VY];
%lambda = diag(Lambda);

G = [A, B; B', -conj(A)];
G = (G + G')/2;
[V, lambda] = eig(G, 'vector');
[lambda, index] = sort(real(lambda), 'descend');
lambda_floor = eps(max(1, max(abs(lambda))))*(2*k);
positive = find(lambda > lambda_floor, k, 'first');
if isempty(positive)
    error('ortho_C_SVQB:Breakdown', ...
        'The input block has no numerically positive C-direction.');
end
lambda = lambda(positive);
V = V(:, index(positive));

VX = V(1 : k, :);
VY = V(k + 1 : end, :);
UX = UX.*d.';
UY = UY.*d.';
U1X = UX*VX + conj_UY.*d.'*VY;
U1Y = UY*VX + conj_UX.*d.'*VY;
sqrt_lambda = sqrt(lambda).';
UX = U1X./sqrt_lambda;
UY = U1Y./sqrt_lambda;
end
