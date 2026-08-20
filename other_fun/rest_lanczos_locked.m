function [Lambda, U, max_res, parout] = rest_lanczos_locked(S, k, Max_time, test, parin)
% This function computes some smallest symplectic eigenvalues of positive definite
% matrices S
% The approach is to used symplectic restarted Lanczos procedure.
% INPUT
%   S: the spd matrix such that H = JS
%   k: the number of smallest eigenvalues of interest
%   p: the addition Lanczos steps performed to improve the accuracy
%   vtil: the initial guess
%   tolcoeff: tolerance in compute the coefficients
%   tolrest: tolerance for stopping restart
%   maxrest: maximal number of times to restart
% OUTPUT
%   Lambda: the k smallest symplectic eigenvalues of S
%   U: symplectic eigenvectors of S

%% default setting
id1 = tic;
[Schol, chS] = chol(S, 'lower');
if chS ~= 0
    error('S is not spd')
end
if ~isfield(parin, 'p');        parin.p = k; end
if ~isfield(parin, 'vtil');     parin.vtil = randn(size(S,1),1); end 
if ~isfield(parin, 'tolc');     parin.tolc = 1e-6; end % tolerance for comupting the coefficients
if ~isfield(parin, 'tolr');     parin.tolr = 1e-7; end % tolr: convergence tolerance
if ~isfield(parin, 'maxr');     parin.maxr = 2*k; end % the maximal number of restartings

p = parin.p;
vtil = parin.vtil;
tolc = parin.tolc;
tolr = parin.tolr;
maxr = parin.maxr;

%% lock parameters
n_locked = 0; % the number of locked
n_active = k; % the number of active(unlocked)
U_locked = []; % locked eigenvectors
V_locked = [];
Lambda_locked = []; % locked eigenvalues
Res = []; % resdiual
nrest = 0;
max_res_locked = 0;

while nrest < maxr
    %% generate Lanczos process of length n_actived + p
    [a, b, V, len] = lanczosing_locked(Schol, n_active + p, vtil, V_locked, tolr, Max_time);
    T = diag(a) + diag(b(2 : end - 1), 1) + diag(b(2 : end - 1),-1);
    [Y, s_new] = eig(T);
    s_new = diag(s_new);
    [s_new, ides] = sort(s_new, 'descend');
    Y = Y(:, ides);
    
    %% the smllest n_active eigenvalues
    n_active = min(n_active, len);
    sT = size(T, 1);
    ei = sqrt(s_new(1 : n_active));
    ei = 1./ei;
    VY = V(:, 1 : sT)*Y(:, 1 : n_active);
    U = horzcat(VY, -Jmul(Schol'\(Schol\VY)));
    U = [-U(:, n_active + 1 : 2*n_active), U(:, 1 : n_active)];
    U1 = -Jmul(U);
    D = sqrt(diag(ei));
    Phi = blkdiag(D, D^-1);
    U = U1*Phi; % the symplectic eigenvectors of M
    In = speye(n_active);
    Qn_active = 1/sqrt(2)*[In, -sqrt(-1)*In;
                           In, sqrt(-1)*In];
    QUQ = test.Qn*U*Qn_active';
    [res_active, res_locked, locked_id] = compute_Res(test.Omega, test.Omega_norm, QUQ, ...
                                        ei, test.n, n_active, tolr);
    
    
    max_res_locked = max(max_res_locked, res_locked);
    max_res = max(res_active,  max_res_locked);

    %% locked
     if ~isempty(locked_id)
        V_new = VY(:, locked_id);
        U_new = U(:, locked_id); % Only half was taken
    
        % % orthonormalize V_new
        % if ~isempty(V_locked)
        %     W_locked = -Jmul(Schol'\(Schol\V_locked));
        %     for i = 1 : size(V_locked, 2)
        %         w_nrm(i) = W_locked(:, i)'*(Schol'\(Schol\W_locked(:, i)));
        %     end 
        %     for i = 1 : length(locked_id)
        %         r = V_locked(:, i)'*(Schol'\(Schol\V_new(:, i)));
        %         V_new(:, i) = V_new(:, i) - r*V_locked(:, i);
        %         nrm = sqrt(V_new(:, i)'*(Schol'\(Schol\V_new(:, i))));
        %         V_new(:, i) = V_new(:, i)/nrm;
        %     end
        % else
        %     for i = 1 : length(locked_id)
        %         nrm = sqrt(V_new(:, i)'*(Schol'\(Schol\V_new(:, i))));
        %         V_new(:, i) = V_new(:, i)/nrm;
        %     end
        % end
        
        V_locked = [V_locked, V_new];
        U_locked = [U_locked, U_new];
        Lambda_locked = [Lambda_locked; ei(locked_id)];
        n_locked = n_locked + length(locked_id);
        n_active = k - n_locked;
     end

    Res = [Res; max_res];

    if n_locked >= k
        break;
    end
    
    %% generate the restarted initinal vector
    if ~isempty(locked_id)
        unlocked_id = setdiff(1 : n_active, locked_id);
        g = coeff(Y(:, unlocked_id), length(unlocked_id), p, tolc);
        vtil = V(:, 1 : sT)*(Y(:, unlocked_id)*g');
    else
        g = coeff(Y(:, 1 : n_active), n_active, p, tolc);
        vtil = V(:, 1 : sT)*(Y(:, 1 : n_active)*g');
    end
    
    if max_res < tolr || toc(id1) > Max_time
        break;
    end
    nrest = nrest + 1;
end
U = [U_locked, U(:, 1 : n_active)];
Lambda = [Lambda_locked; ei(1:n_active)];
[Lambda, id] = sort(Lambda);
U = U(:, id);

if nrest == maxr
    warning('Restarting seems to fail')
    mess = 'Restarting seems to fail';
else
    mess = 'Restarting successfully';
end
parout = struct;
parout.mess = mess;
parout.nrest = nrest;
end