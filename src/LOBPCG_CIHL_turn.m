function [RX, RY, ZX, ZY, Theta_pos, Res, lock_ID, time_mat] = LOBPCG_CIHL_turn(A,...
                                            B, conj_A, conj_B, X0, Y0, T, par)
%% Compute the l_pos smallest eigenpairs of the pair (Omega, Cn)
% under the Cn-inner product using the improved HL trick
id4 = tic;
tol = par.tol; % convergent tolerance
Max_iter = par.Max_iter;  % maximal iteration
Max_time = par.Max_time; % maximal time
l_pos = par.l_pos; % the number of desired positive eigenvalues
k_pos = par.k_pos; % more than l_pos
n = par.n;
criterion = par.criterion;
if ~isfield(par, 'Omega_norm')
    Rand = randn(n, 10);
    OmegaRand = [A*Rand + B*Rand; conj_B*Rand + conj_A*Rand];
    Omega_norm = norm(OmegaRand, 'fro')/norm(Rand, 'fro');
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

PX = []; PY = []; Res = zeros(Max_iter, 1);

arg = struct; 
arg.sort = 0;
arg.stable = 0;
arg.reorth = 0;

%% Cn-Orthonormalize X
[UX, UY] = ortho_C_SVQB(X0, Y0);
[UX, UY] = ortho_C_SVQB(UX, UY);

%% Rayleigh-Ritz procedure
OmegaUX = A*UX + B*UY;
OmegaUY = conj_B*UX + conj_A*UY;
A_k = UX'*OmegaUX + UY'*OmegaUY;
B_k = UX'*conj(OmegaUY) + UY'*conj(OmegaUX);
A_k = (A_k + A_k')/2;
B_k = (B_k + B_k.')/2;
[VX, VY, Theta_pos] = BSE_complex(A_k, B_k);
index = 1 : k_pos;
Theta_pos = Theta_pos(index, index);
theta_pos = diag(Theta_pos);
ZX = UX*VX + conj(UY)*VY;
ZY = UY*VX + conj(UX)*VY;
Z = [ZX; ZY];
Exit = false;
star_check = 0;
np = 0;
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
    Res(i) = res;
   
    if res < 1e-10
        star_check = 1;
    end

    if star_check && i > 5
        if criterion == 4 && i <= 10
         Exit = false;
        else
         Exit = CheckExit(log10(Res(1:i)), i, criterion);
        end
    end
    
    if res < tol || Exit == true || toc(id4) > Max_time
        RX = RX(:, 1 : k_pos);
        RY = RY(:, 1 : k_pos);
        ZX = ZX(:, 1 : k_pos);
        ZY = ZY(:, 1 : k_pos);
        Theta_pos = Theta_pos(1 : k_pos, 1 : k_pos);
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

    %% deflate P  
    if i > 1
        activeID_pos = activeID_pos(activeID_pos <= size(PX, 2));
        PX = PX(:, activeID_pos);
        PY = PY(:, activeID_pos);
        np = size(activeID_pos, 2);
    end

    ZPX = [ZX, PX];
    ZPY = [ZY, PY];

    %% enhance orthogonalization
    if test_Cortho(ZPX, ZPY) > min(1e-8, 0.1*res)
        arg.stable = 0;
        [ZPX, ZPY, CX, CY] = ortho_C_CGS(ZPX, ZPY, arg);
        CX = CX(1:k_pos, 1:k_pos);
        CY = CY(1:k_pos, 1:k_pos);
        oldOmegaZX = OmegaZX;
        OmegaZX = oldOmegaZX*CX + conj(OmegaZY)*CY;
        OmegaZY = OmegaZY*CX + conj(oldOmegaZX)*CY;
    end

    %% protection
    %if test_ortho(ZPX, ZPY) > 1e-8
    %    fprintf("Trigger protection at %dth step\n", i);
    %    arg.stable = 1;
    %    [ZPX, ZPY] = ortho_C_CGS(ZPX, ZPY, arg);
    %end

    %% W is orthonormalize to [Z P]
    p1 = size(ZPX, 2);
    p2 = size(WX, 2 );
    if(p1 + p2 > n)
        ZPX = ZPX(:, 1 : n - p2);
        ZPY = ZPY(:, 1 : n - p2);
        np = n - p2 - k_pos;
    end
    [WX, WY] = ortho_together_C(WX, WY, ZPX, ZPY);

    %% Project
    UX = [ZPX, WX];
    UY = [ZPY, WY];

    [Theta, VX, VY, time_mat] = RayleighRitz_C(UX, UY, OmegaZX, OmegaZY, A, B, conj_A, ...
                                        conj_B, k_pos, np, time_mat, 'aug');
    Theta_pos = Theta(1 : k_pos, 1 : k_pos);
    theta_pos = diag(Theta_pos);

    %% Update P, Z
    [ZX, ZY, PX, PY] = update_CIHL(UX, UY, VX, VY, k_pos);
    Z = [ZX; ZY];
end
Res = Res(1:i);
RX = RX(:, 1 : k_pos);
RY = RY(:, 1 : k_pos);
ZX = ZX(:, 1 : k_pos);
ZY = ZY(:, 1 : k_pos);
Theta_pos = Theta_pos(1 : k_pos, 1 : k_pos);
