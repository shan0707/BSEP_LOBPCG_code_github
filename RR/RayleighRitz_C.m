function [Theta, VX, VY, time_mat] = RayleighRitz_C(UX, UY, OmegaZX, OmegaZY,...
                                A, B, conj_A, conj_B, k_pos, np, time_mat, arg)
%% Rayleigh-Ritz with Cn-inner product
%output the k_pos smallest positive eigenpairs of (U'OmegaU, Ck)
%UOmegaU = [A_k, B_k; conj(B_k), conj(A_k)]

[A_k, B_k, time_mat] = generate_UOmegaU(UX, UY, OmegaZX, OmegaZY, A, B, ...
                                        conj_A, conj_B, k_pos, np, time_mat);


if nargin > 11
    [VX, VY, Theta] = BSE_complex(A_k, B_k); % output all positive eigenpairs
else 
    [VX, VY, Theta] = BSE_complex(A_k, B_k, k_pos); % output k_pos positive eigenpairs
end