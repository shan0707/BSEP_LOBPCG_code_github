function [UX, UY, time_U] = ortho_Omega_BSE(UX, UY, Omega, n)
% Omega-orthonormalize [UX, conj(UY); UY, conj(UX)] by solving a BSEP
id1 = tic;
OmegaU = Omega*[UX; UY];
time_U = toc(id1);

OmegaUX = OmegaU(1 : n, :);
OmegaUY = OmegaU(n + 1 : end, :);

A_k = UX'*OmegaUX + UY'*OmegaUY;
B_k = UX'*conj(OmegaUY) + UY'*conj(OmegaUX);

d = abs(diag(A_k)).^(-1/2);
A_k = d.*A_k.*d.';
B_k = d.*B_k.*d.';
A_k = (A_k + A_k')/2;
B_k = (B_k + B_k.')/2;

[VX, VY, Lambda] = BSE_complex(A_k, B_k);

UX = UX.*d.';
UY = UY.*d.';
U1X = UX*VX + conj(UY)*VY;
U1Y = UY*VX + conj(UX)*VY;
sqrt_lambda = diag(Lambda).^(-1/2);
UX = U1X.*sqrt_lambda.';
UY = U1Y.*sqrt_lambda.';