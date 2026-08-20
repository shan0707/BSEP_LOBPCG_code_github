addpath(genpath(pwd));
clear;
rng(0);
delete('./Data/Sparse*.txt');

MatNames = ["bcsstk21"; "fv1"; "crystm03"; "apache1"; "shallow_water2"];
Dimension = [3600, 9604, 24696, 80800, 81920];
num_test = 5; % the number of test matrices
num_alg = 6;
list_k = [18, 48, 100, 100, 100]; % the number of computed eigenvalues
Max_time = 3600;
Max_iter = 200;

time = zeros(num_alg, num_test);
rest = zeros(num_alg, num_test); % record the residual;
eig_cell = cell(num_alg, num_test);

Res_CIHL = zeros(num_test, Max_iter);
Res_CIHL_OmegaIHL = zeros(num_test, Max_iter);
Iter = zeros(1, num_test);

for i = 1 : num_test
    %% generate test matrix symplectic Gauss transformation M
    disp(MatNames(i));
    FileName = MatNames(i) + '.mat';
    load(FileName);
    M = Problem.A;
    k = list_k(i);
    n = Dimension(i)/2;
    In = speye(n);
    Qn = 1/sqrt(2)*[In, -sqrt(-1)*In;
                    In, sqrt(-1)*In];
    Ik = speye(k);
    Qk = [Ik, -sqrt(-1)*Ik;
          Ik, sqrt(-1)*Ik];

    H_ham = sparse([M(n + 1 : 2*n, :); -M(1 : n, :)]); % H_ham: hamiltonian matrix
    H = Qn*(-sqrt(-1)*H_ham)*Qn'; % H: BSE matrix
    Omega = Qn*M*Qn';
    A = Omega(1 : n, 1 : n);
    B = Omega(1 : n, n + 1 : 2*n);
    
    % parameters for test
    test = struct;
    test.Qk = Qk;
    test.Omega = Omega;
    test.Qn = Qn;
    test.n = n;
    Omega_norm = normest(M, 1e-4);
    test.Omega_norm = Omega_norm;
   
    %% restart Lanczos
    par = struct;
    par.p = max(k, 50 - k); % parameter for restarting
    par.tolr = 1e-14 ; % tolerance for restarting
    par.maxr = 1000;  % the maximal number of restartings
    par.tolc = 1e-12; % tolerance for comupting the coefficients
    out_splLanz = cell(1, num_test); % restarted symplectic lanczos

    disp('Warming up...');
    rest_lanczos_locked(M, k, 50, test, par);
    disp('Warm-up finished.');

    ID1 = tic;
    [eig_cell{1, i}, ~, res_lanczos, ~] = rest_lanczos_locked(M, k, Max_time, test, par);
    time(1, i) = toc(ID1);
    
    %% Riemannian optimizing
    opts = struct;
    opts.record = 0;
    opts.mxitr = 1000000;
    opts.xtol = 1e-11;
    opts.ftol = 1e-11;
    opts.gtol = 1e-10;
    opts.maxtau = 1;
    out_spopt = cell(1, num_test);
    opts.tol = 1e-14;

    X0 = zeros(2*n, 2*k);
    X0(1 : k, 1 : k) = eye(k);
    X0(n + 1 : n + k, k + 1 : end) = eye(k);

    disp('Warming up...');
    spopt(X0, @eigvalcost, opts, 100, test, M);
    disp('Warm-up finished.');

    ID2 = tic;
    [~, res_Rieman, eig_cell{2, i}, ~]= spopt(X0, @eigvalcost, opts, Max_time, test, M);
    time(2, i) = toc(ID2);

    %% LOBPCG
    lob = struct;
    lob.tol = 1e-14 ; % tolerance
    lob.Max_time = Max_time; % maximal time
    lob.l_pos = k; % the number of desired eigenvalues
    lob.k_pos = max(ceil(1.5*k), k + 5);
    lob.n = n; % the dimension of the matrix A
    lob.Omega_norm = Omega_norm; % the norm of Omega
    lob.pre = 3; % the way of computing preconditioner

    X0 = randn(n, lob.k_pos);
    Y0 = randn(n, lob.k_pos);
    T = blkdiag(A, conj(A)); % preconditioner
    T = sparse(T);

    %% LOBPCG-CIHL
    lob.Max_iter = 2; % maximal iteration
    disp('Warming up...');
    LOBPCG_CIHL(Omega, X0, Y0, T, lob);
    disp('Warm-up finished.');

    ID3 = tic;
    lob.Max_iter = 200;  % maximal iteration
    [~, ~, ~, ~, Theta, Res, ~, time_mat] = LOBPCG_CIHL(Omega, X0, Y0, T, lob);
    time(3, i) = toc(ID3)
    CIHL_mat = time_mat/time(3, i)
    res_CIHL = Res(end)
    eig_cell{3, i} = diag(Theta(1 : k, 1 : k));
    iter_CIHL = length(Res)
    Res_CIHL(i, 1 : iter_CIHL) = Res;
        

    %% LOBPCG-CIHL-OmegaIHL
    lob.reorth_Omega = 1; % the number of reorthogonalization
    lob.criterion = 4; % turning algorithm criterion
    lob.Max_iter = 2;
    disp('Warming up...');
    LOBPCG_CIHL_OmegaIHL(Omega, X0, Y0, T, lob);
    disp('Warm-up finished.');

    ID4 = tic;
    lob.Max_iter = 200;
    [~, ~, ~, ~, Theta, Res, ~, time_mat, arg] = LOBPCG_CIHL_OmegaIHL(Omega, X0, Y0, T, lob);
    time(4, i) = toc(ID4)
    CIHL_OmegaIHL_mat = time_mat/time(4, i)
    res_CIHL_OmegaIHL = Res(end)
    eig_cell{4, i} = diag(Theta(1 : k, 1 : k));
    iter_CIHL_OmegaIHL = length(Res)
    Res_CIHL_OmegaIHL(i, 1 : iter_CIHL_OmegaIHL) = Res;

   
    nk = [Dimension(i), k];
    Iter_Res_Time = [0, 0, iter_CIHL, iter_CIHL_OmegaIHL;
                     res_lanczos, res_Rieman, res_CIHL, res_CIHL_OmegaIHL;
                     time(1, i), time(2, i), time(3, i), time(4, i)];
   
    dlmwrite('./Data/Sparse_iter_res_time.txt', nk, '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite('./Data/Sparse_iter_res_time.txt', Iter_Res_Time, '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite('./Data/Sparse_rest.txt', Res_CIHL(i, :), '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite('./Data/Sparse_rest.txt', Res_CIHL_OmegaIHL(i, :), '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite('./Data/Sparse_iter.txt', iter_CIHL, '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite('./Data/Sparse_iter.txt', iter_CIHL_OmegaIHL, '-append', 'delimiter', ',', 'precision', 4);
end

%% ============= functions ==========================
function [F,G] = eigvalcost(X,M)
MX = M*X;
% F, G represent object function and gradient
F = X(:)'*MX(:);
G = 2*MX;
end