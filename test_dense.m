addpath(genpath(pwd));
clear;
rng(0);
delete('./Data/dense*.txt')

num_test = 5; % the number of test matrices
num_alg = 4;
L_pos = [3, 12, 23, 50, 50]; % the number of computed eigenvalues
Max_time = 1800;
Max_iter = 200;

time = zeros(num_alg, num_test);
rest = zeros(num_alg, num_test); % record the residual;
eig_cell = cell(num_alg, num_test);

Res_Davidson = zeros(num_test, Max_iter);
Res_SD = zeros(num_test, Max_iter);
Iteration = [];

for i = 1 : num_test
    %% generate test matrix
    [A, B, n] = generate_BSE(i);
    disp(n);
    Omega = [A, B; conj(B), conj(A)];
    In = speye(n);
    Qn = 1/sqrt(2)*[In, -sqrt(-1)*In;
                    In, sqrt(-1)*In];
    k = L_pos(i);
    Ik = speye(k);
    Qk = [Ik, -sqrt(-1)*Ik;
          Ik, sqrt(-1)*Ik];
    M = Qn'*Omega*Qn;

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
    par.p = min(n, par.p);
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
    opts.mxitr = 2000000;
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
    
    %%
    par = struct;
    par.tol = 1e-14 ; % tolerance
    par.Max_iter = Max_iter;  % maximal iteration
    par.Max_time = Max_time; % maximal time
    par.l_pos = k; % the number of desired eigenvalues
    par.k_pos = max(ceil(1.5*k), k + 5);
    par.n = n; % the dimension of the matrix A
    par.Omega_norm = Omega_norm; % the norm of Omega
    par.reorth_Omega = 1; % the number of reorthogonalization

    X0 = randn(n, par.k_pos);
    Y0 = randn(n, par.k_pos);
    a = diag(A);
    T = blkdiag(diag(a), diag(conj(a)));
    T = sparse(T);

    %% Davidson
    name = "Davidson";
    disp(name);
    par.Max_iter = 2;
    disp('Warming up...');
    Davidson(Omega, X0, Y0, par);
    disp('Warm-up finished.');

    ID3 = tic;
    par.Max_iter = 40000;
    [~, ~, ~, ~, Theta3, Res3, ~, ~] = Davidson(Omega, X0, Y0, par);
    time(3, i) = toc(ID3)
    iter_Davidson = length(Res3);
    res_Davidson = Res3(end);
    eig_cell{3, i} = diag(Theta3(1 : k, 1 : k));
    Res_Davidson(i, 1 : iter_Davidson) = Res3;
    
    %% Steepest_Descent
    name = "Steepest_Descent";
    disp(name);n
    par.Max_iter = 2;
    disp('Warming up...');
    Steepest_Descent(Omega, X0, Y0, T, par);
    disp('Warm-up finished.');
    
    ID4 = tic;
    par.Max_iter = 20000;
    [~, ~, ~, ~, Theta4, Res4, ~, ~] = Steepest_Descent(Omega, X0, Y0, T, par);
    time(4, i) = toc(ID4)
    iter_SD = length(Res4);
    res_SD = Res4(end);
    eig_cell{4, i} = diag(Theta4(1 : k, 1 : k));
    Res_SD(i, 1 : iter_SD) = Res4;


    nk = [n, k];
    Iter_Res_Time = [0, 0,   iter_Davidson, iter_SD;
                     res_lanczos, res_Rieman, res_Davidson, res_SD;
                     time(1, i), time(2, i), time(3, i), time(4, i)];
   
    dlmwrite('./Data/dense_iter_res_time.txt', nk, '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite('./Data/dense_iter_res_time.txt', Iter_Res_Time, '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite('./Data/dense_rest.txt', Res_Davidson(i, :), '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite('./Data/dense_rest.txt', Res_SD(i, :), '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite('./Data/dense_iter.txt', iter_Davidson, '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite('./Data/dense_iter.txt', iter_SD, '-append', 'delimiter', ',', 'precision', 4);
end

hold off;

%% ============= functions ==========================
function [F,G] = eigvalcost(X,M)
MX = M*X;
% F, G represent object function and gradient
F = X(:)'*MX(:);
G = 2*MX;
end