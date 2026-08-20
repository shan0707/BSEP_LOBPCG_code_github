addpath(genpath(pwd));
clear;
rng(1);
delete('./Data/ex2*.txt');

k = 5; % the number of desired eigenvalues 
list_n = [2000, 5000];
ln = length(list_n); % the number of test matrices
Max_time = 1800;
num_alg = 3;

time = zeros(num_alg, ln);
rest = zeros(num_alg, ln); % record the residual;
eig_cell = cell(num_alg, ln); 

Ik = speye(k);
Qk = [Ik, -sqrt(-1)*Ik;
      Ik, sqrt(-1)*Ik];

% parameters for Riemannian optimization 
opts = struct; 
opts.record = 0;
opts.mxitr = 1000000;
opts.xtol = 1e-11;
opts.ftol = 1e-11;
opts.gtol = 1e-10;
opts.maxtau = 1;
out_spopt = cell(1, ln);
opts.tol = 1e-14;

% parameters for restarted symplect is lanczos
par = struct;
par.p = max(k, 50 - k); % parameter for restarting
par.tolr = 1e-14; % tolerance for restarting
par.maxr = 100000;  % the maximal number of restartings
par.tolc = 1e-12; % tolerance for comupting the coefficients
out_splLanz = cell(1, ln); % restarted symplectic lanczos

% parameters for LOBPCG
lob = struct;
lob.tol = 1e-14 ; % tolerance
lob.Max_iter = 200;  % maximal iteration
lob.Max_time = Max_time; % maximal time

for i = 1 : ln
    %% generate test matrix symplectic Gauss transformation M
    n = list_n(i);
    disp(n);
    [M, ~] = generate_ex2(n);

    In = speye(n);
    Qn = 1/sqrt(2)*[In, -sqrt(-1)*In;
                    In, sqrt(-1)*In];
 
    H_ham = sparse([M(n+1:2*n,:); -M(1:n,:)]); % H_ham: hamiltonian matrix
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
    Rand = randn(2*n, k);
    Omega_norm = norm(Omega*Rand, 'fro')/norm(Rand, 'fro');
    test.Omega_norm = Omega_norm;
    
    %% restart Lanczos
    disp('Warming up...');
    rest_lanczos_locked(M, k, 50, test, par);
    disp('Warm-up finished.');

    ID1 = tic;
    [eig_cell{1, i}, ~, max_res, out_splLanz{i}] = ...
                            rest_lanczos_locked(M, k, Max_time, test, par);
    time(1, i) = toc(ID1);
    rest(1, i) = max_res;
    
    %% Riemannian optimizing
    X0 = zeros(2*n, 2*k); 
    X0(1 : k, 1 : k) = eye(k);
    X0(n + 1 : n + k, k + 1 : end) = eye(k);

    disp('Warming up...');
    spopt(X0, @eigvalcost, opts, 50, test, M);
    disp('Warm-up finished.');

    ID2 = tic;
    [Xmin, max_res, eig_cell{2, i}, out_spopt{i}] = ...
                           spopt(X0, @eigvalcost, opts, Max_time, test, M);
    time(2, i) = toc(ID2);
    rest(2, i) = max_res;

    %% LOBPCG-CIHL-OmegaIHL
    name = "LOBPCG-CIHL-OmegaIHL";
    disp(name);
    lob.n = n;
    lob.l_pos = k;
    lob.k_pos = max(ceil(1.5*k), k + 5);
    lob.Omega_norm = Omega_norm;
    
    X0 = randn(n, lob.k_pos);
    Y0 = randn(n, lob.k_pos); 
    T = [diag(diag(A)), diag(diag(B)); 
         diag(conj(diag(B))), diag(diag(conj(A)))]; %preconditioner
    T = sparse(T);
    
    lob.reorth_Omega = 1;
    lob.Max_iter = 2;
    lob.criterion = 4;
    disp('Warming up...');
    LOBPCG_CIHL_OmegaIHL(Omega, X0, Y0, T, lob);
    disp('Warm-up finished.');

    ID3 = tic;
    lob.Max_iter = 200;
    [~, ~, ~, ~, Theta3, Res3, ~, time_mat, arg] = LOBPCG_CIHL_OmegaIHL(Omega, X0, Y0, T, lob);
    time(3, i) = toc(ID3)
    CIHL_OmegaIHL_mat = time_mat/time(3, i)
    rest(3, i) = Res3(end)
    eig_cell{3, i} = diag(Theta3(1 : k, 1 : k));
end

dlmwrite('./Data/ex2_time.txt', time, '-append', 'delimiter', ',', 'precision', 4);
dlmwrite('./Data/ex2_rest.txt', rest, '-append', 'delimiter', ',', 'precision', 4);
dlmwrite('./Data/ex2_res.txt', Res3, '-append', 'delimiter', ',', 'precision', 4);
dlmwrite('./Data/ex2_eig.txt', eig_cell, '-append', 'delimiter', ',', 'precision', 15);



%% ============= functions ==========================
function [F,G] = eigvalcost(X,M)
MX = M*X;
% F, G represent object function and gradient
F = X(:)'*MX(:);
G = 2*MX;
end