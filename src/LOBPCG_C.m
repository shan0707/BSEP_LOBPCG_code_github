function [RX, RY, ZX, ZY, Theta_pos, Res, lock_ID, time_mat] = LOBPCG_C(Omega, X0, Y0, T, par)
%% Compute the l_pos smallest eigenpairs of the pair (Omega, Cn) under the Cn-inner product
id1 = tic;
tol = par.tol; % convergent tolerance
Max_iter = par.Max_iter;  % maximal iteration
Max_time = par.Max_time; % maximal time
l_pos = par.l_pos; % the number of desired positive eigenvalues
k_pos = par.k_pos; % more than l_pos
n = par.n;
if ~isfield(par, 'Omega_norm')
    Rand = randn(2*n, 10);
    Omega_norm = norm(Omega*Rand, 'fro')/norm(Rand, 'fro');
else
    Omega_norm = par.Omega_norm;
end
if ~isfield(par, 'pre')
    pre = 1;
else
    pre = par.pre;
end
if pre ~= 1
    [TL, TU] = preconditioner_LU(T, pre);
end

A = Omega(1 : n, 1 : n);
B = Omega(1 : n, n + 1 : end);
conj_B = Omega(n + 1 : 2*n, 1 : n);
conj_A = Omega(n + 1 : 2*n, n + 1 : 2*n);
PX = []; PY = []; Res = [];

%% Cn-Orthonormalize X
arg = struct; 
arg.sort = 0;
arg.stable = 0;
arg.reorth = 1;
[UX, UY] = ortho_C_CGS(X0, Y0, arg);

%% Rayleigh-Ritz procedure 
OmegaUX = A*UX + B*UY;
OmegaUY = conj_B*UX + conj_A*UY;
A_k = UX'*OmegaUX + UY'*OmegaUY;
B_k = UX'*conj(OmegaUY) + UY'*conj(OmegaUX);
A_k = (A_k + A_k')/2; B_k = (B_k + B_k.')/2;
[VX, VY, Theta_pos] = BSE_complex(A_k, B_k);
index = 1 : k_pos;
Theta_pos = Theta_pos(index, index);
ZX = UX*VX + conj(UY)*VY;
ZY = UY*VX + conj(UX)*VY;
Z = [ZX; ZY];
np = 0;
time_mat = 0;

for i = 1 : Max_iter
    id2 = tic;
    OmegaZX = A*ZX + B*ZY;
    OmegaZY = conj_B*ZX + conj_A*ZY;
    time_mat = time_mat + toc(id2);
    RX = OmegaZX - ZX*Theta_pos;
    RY = OmegaZY + ZY*Theta_pos;
    R = [RX; RY];

    %% Deflation
    [res, lockID_pos, lockID_neg, ~] = deflation(R, Theta_pos, Z, ...
                                            Omega_norm, k_pos, l_pos, tol);
    lock_ID = [lockID_pos, lockID_neg];
    Res = [Res; res];
    
    if res < tol || toc(id1) > Max_time
        RX = RX(:, 1 : l_pos);
        RY = RY(:, 1 : l_pos);
        ZX = ZX(:, 1 : l_pos);
        ZY = ZY(:, 1 : l_pos);
        Theta_pos = Theta_pos(1 : l_pos, 1 : l_pos);
        break;
    end
    activeID_pos = setdiff((1 : k_pos), lockID_pos);

    %% preconditioner
    if pre == 1
        W = T\R(:, activeID_pos);
    else
        W = preconditioner_solve(TL, TU, R(:, activeID_pos), pre);
    end

    WX = W(1 : n, :);
    WY = W(n + 1 : end, :);

    if i > 1
        activeID_pos = activeID_pos(activeID_pos <= size(PX, 2));
        PX = PX(:, activeID_pos);
        PY = PY(:, activeID_pos);
        np = size(activeID_pos, 2);
    end

    %% Project
    UX = [ZX, WX, PX];
    UY = [ZY, WY, PY];
    
    arg.stable = 0;
    [UX, UY] = ortho_C_CGS(UX, UY, arg);

    %% protection
    if test_Cortho(UX, UY) > 1e-8
        fprintf("Trigger protection at %dth step\n", i);
        arg.stable = 1;
        [UX, UY] = ortho_C_CGS(UX, UY, arg);
    end

    [Theta_pos, VX, VY, time_mat] = RayleighRitz_C(UX, UY, OmegaZX, OmegaZY, ...
                                A, B, conj_A, conj_B, k_pos, np, time_mat);

    %% Update P, Z
    [ZX, ZY, PX, PY] = update_C(UX, UY, VX, VY, k_pos);

    Z = [ZX; ZY];
end
RX = RX(:, 1 : l_pos);
RY = RY(:, 1 : l_pos);
ZX = ZX(:, 1 : l_pos);
ZY = ZY(:, 1 : l_pos);
Theta_pos = Theta_pos(1 : l_pos, 1 : l_pos);