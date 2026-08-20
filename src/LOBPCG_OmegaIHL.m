function [RX, RY, ZX, ZY, Theta_pos, Res, lock_ID, time_mat] = LOBPCG_OmegaIHL(Omega, X0, Y0, T, par)
%% Compute the l_pos smallest eigenpairs of the pair (Omega, Cn) under 
% the Omega-inner product using the improved HL trick
id2 = tic;
tol = par.tol; % convergent tolerance
Max_iter = par.Max_iter;  % maximal iteration
Max_time = par.Max_time; % maximal time
l_pos = par.l_pos; % the number of desired positive eigenvalues
k_pos = par.k_pos; % more than l_pos
n = par.n;
reorth_Omega = par.reorth_Omega;
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
conj_A = Omega(n + 1 : 2*n,  n + 1 : 2*n);
conj_B = Omega(n + 1 : 2*n, 1 : n);
PX = []; 
PY = [];
OmegaPX = [];
OmegaPY = [];
Res = [];

%% Omega-Orthonormalize X
[UX, UY] = ortho_Omega_CGS(X0, Y0, Omega, reorth_Omega);

%% Rayleigh-Ritz procedure
[Theta_pos, VX, VY] = RayleighRitz_Omega(UX, UY);
theta_pos = 1./diag(Theta_pos);
Theta_pos = diag(theta_pos);
ZX = UX*VX + conj(UY)*VY;
ZY = UY*VX + conj(UX)*VY;
Z = [ZX; ZY];
time_mat = 0;

for i = 1 : Max_iter
    id1 = tic;
    OmegaZX = A*ZX + B*ZY;
    OmegaZY = conj_B*ZX + conj_A*ZY;
    time_mat = time_mat + toc(id1);
    RX = OmegaZX - ZX.*theta_pos.';
    RY = OmegaZY + ZY.*theta_pos.';
    R = [RX; RY];
    
    %% Deflation
    [res, lockID_pos, lockID_neg, ~] = deflation(R, Theta_pos, Z, ...
                                            Omega_norm, k_pos, l_pos, tol);
    lock_ID = [lockID_pos, lockID_neg];
    Res = [Res; res];

    if res < tol || toc(id2) > Max_time
        RX = RX(:, 1 : l_pos);
        RY = RY(:, 1 : l_pos);
        ZX = ZX(:, 1 : l_pos);
        ZY = ZY(:, 1 : l_pos);
        Theta_pos = Theta_pos(1 : l_pos, 1 : l_pos);
        break;
    end
    activeID_pos = setdiff((1 : k_pos), lockID_pos); % unlocked index
   
    %% preconditioner
    if pre == 1
        W = T\R(:, activeID_pos);
    else
        W = preconditioner_solve(TL, TU, R(:, activeID_pos), pre);
    end
    WX = W(1 : n, :);
    WY = W(n + 1 : end, :);

    ZPX = [ZX, PX];
    ZPY = [ZY, PY];

    %% W is orthonormalize to [Z P]
    p1 = size(ZPX, 2);
    p2 = size(WX, 2);
    if(p1 + p2 > n)
        ZPX = ZPX(:, 1 : n - p2);
        ZPY = ZPY(:, 1 : n - p2);
    end
    if i > 1
        id1 = tic;
        OmegaPX = A*PX + B*PY;
        OmegaPY = conj_B*PX + conj_A*PY;
        time_mat = time_mat + toc(id1);
    end
    OmegaZPX = [OmegaZX, OmegaPX];
    OmegaZPY = [OmegaZY, OmegaPY];
    [WX, WY, time_mat] = ortho_together_Omega(WX, WY, ZPX, ZPY, ...
                                        OmegaZPX, OmegaZPY, Omega, time_mat);
    UX = [ZPX, WX];
    UY = [ZPY, WY];
   
    %% Project
    [Theta, VX, VY] = RayleighRitz_Omega(UX, UY);
    Theta_pos = Theta(1 : k_pos, 1 : k_pos);
    theta_pos = 1./diag(Theta_pos);
    Theta_pos = diag(theta_pos);
    
    %% Update P, Z
    [ZX, ZY, PX, PY] = update_OmegaIHL(UX, UY, VX, VY, k_pos);
    Z = [ZX; ZY];
end
RX = RX(:, 1 : l_pos);
RY = RY(:, 1 : l_pos);
ZX = ZX(:, 1 : l_pos);
ZY = ZY(:, 1 : l_pos);
Theta_pos = Theta_pos(1 : l_pos, 1 : l_pos);