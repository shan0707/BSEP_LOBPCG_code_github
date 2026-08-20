function [UX, UY, CX, CY] = ortho_C_CGS(UX, UY, arg)
%% Cn-orthonormalize [UX, conj(UY); UY, conj(UX)] use the CGS
% Arguments:
% UX, UY  - Input basis matrices to be orthonormalized.
% arg.sort = 1; % 1-sort: performs pre-conditioning by sorting columns
% arg.stable = 1; % 1: activates the high-stability mode using explicit projection operators.
% arg.reorth = 1; % Number of re-orthogonalization steps (1 or 2).
% [UX, conj(UY); UY, conj(UX)] = [UX_old, conj(UY_old); UY_old, conj(UX_old)]*[CX, conj(CY); CY, conj(CX)];

[n, p] = size(UX);
steps = 1 + arg.reorth;
CX = eye(p, 'like', UX);
CY = zeros(p, 'like', UX);

%% Pre-conditioning: Sort columns to enhance numerical stability
if arg.sort == 1
    D = real(sum(abs(UX).^2, 1) - sum(abs(UY).^2, 1));
    [~, index] = sort(D, 'descend');
    UX = UX(:, index);
    UY = UY(:, index);
    CX = CX(:, index);
    CY = CY(:, index);
end

%% Cn-orthogonalization via CGS
for i = 1 : p
    if i == 1
        [UX(:, 1), UY(:, 1), CX(:, 1), CY(:, 1)] = ...
            C_orthonormalize(UX(:, 1), UY(:, 1), CX(:, 1), CY(:, 1));
    else
        ui = [UX(:, i); UY(:, i)];
        alpha = CX(:, i);
        beta = CY(:, i);
        index = 1 : i - 1;
        UX_id = UX(:, index);
        UY_id = UY(:, index);

        conjUX_id = conj(UX_id);
        conjUY_id = conj(UY_id);

        if arg.stable == 0
            % Fast Mode: Inline block multiplications
            for k = 1 : steps
                a = ui(1 : n);
                b = ui(n + 1 : end);

                % Step 1: equivalent to CUH_multi
                c_tmp = UX_id'*a - UY_id'*b;
                d_tmp = conjUY_id'*a - conjUX_id'*b;

                % Step 2: equivalent to UC_multi
                c_out = UX_id*c_tmp - conjUY_id*d_tmp;
                d_out = UY_id*c_tmp - conjUX_id*d_tmp;

                ui = ui - [c_out; d_out];
                alpha = alpha - CX(:, index)*c_tmp + conj(CY(:, index))*d_tmp;
                beta = beta - CY(:, index)*c_tmp + conj(CX(:, index))*d_tmp;
            end
        else
            % High-Stability Mode (Explicit matrices)
            CnU = [UX_id,  conjUY_id; -UY_id, -conjUX_id];
            UCk = [UX_id, -conjUY_id;  UY_id, -conjUX_id];
            for k = 1 : steps
                coefficient = CnU'*ui;
                c_tmp = coefficient(1:i-1);
                d_tmp = coefficient(i:end);
                ui = ui - UCk*coefficient;
                alpha = alpha - CX(:, index)*c_tmp + conj(CY(:, index))*d_tmp;
                beta = beta - CY(:, index)*c_tmp + conj(CX(:, index))*d_tmp;
            end
        end

        uix = ui(1 : n, :);
        uiy = ui(n + 1 : end, :);
        [UX(:, i), UY(:, i), CX(:, i), CY(:, i)] = ...
            C_orthonormalize(uix, uiy, alpha, beta);
    end
    % the number of orthogonal columns exceeds the numer of rows
    if i == n
        UX = UX(:, 1 : n);
        UY = UY(:, 1 : n);
        CX = CX(:, 1 : n);
        CY = CY(:, 1 : n);
        break;
    end
end



function [ux, uy, alpha_out, beta_out] = C_orthonormalize(x, y, alpha, beta)
% Cn-orthonormalize vectors and maintain the required symmetry structure.
% If the inner product is negative, swap and conjugate to preserve
% the deflating subspace properties.
inner_prod = real(x'*x - y'*y);
norm_factor = sqrt(abs(inner_prod));
ux = x/norm_factor;
uy = y/norm_factor;
alpha_out = alpha/norm_factor;
beta_out = beta/norm_factor;

if inner_prod < 0
    temp = ux;
    ux = conj(uy);
    uy = conj(temp);
    temp = alpha_out;
    alpha_out = conj(beta_out);
    beta_out = conj(temp);
end
