addpath(genpath(pwd));
clear;
rng(1);
delete('./Data/LOBPCG_iter_res_time.txt');
delete('./Data/Iter_case*.txt')
delete('./Data/Res_case*.txt')

num_test = 5; % the number of test matrices
num_alg = 6;
L_pos = [3, 12, 23, 50, 50]; % the number of computed eigenvalues

% parameters for LOBPCG
par = struct;
par.tol = 1e-14 ; % tolerance
par.Max_iter = 200;  % maximal iteration
par.Max_time = 1800; % maximal time

time = zeros(num_alg, num_test);
rest = zeros(num_alg, num_test); % record the residual;
eig_cell = cell(num_alg, num_test);

Res_C = zeros(num_test, par.Max_iter);
Res_OmegaIHL = zeros(num_test, par.Max_iter);
Res_CIHL = zeros(num_test, par.Max_iter);
Res_CIHL_OmegaIHL = zeros(num_test, par.Max_iter);
Iteration = [];

for i = 1 : num_test
    %% generate test matrix
    [A, B, n] = generate_BSE(i);
    Omega = [A, B; conj(B), conj(A)];

    disp(n);
    par.n = n;
    par.l_pos = L_pos(i);
    par.k_pos = max(ceil(1.5*par.l_pos), par.l_pos + 5);
  
    %% setting
    X0 = randn(n, par.k_pos);
    Y0 = randn(n, par.k_pos);
    a = diag(A);
    T = spdiags([a; conj(a)], 0, 2*n, 2*n); % preconditioner
    Rand = randn(2*n, 10);
    par.Omega_norm = norm(Omega*Rand, 'fro')/norm(Rand, 'fro');
 
    %% LOBPCG-C
     name = "LOBPCG-C";
     disp(name);
     par.Max_iter = 2;
     disp('Warming up...');
     LOBPCG_C(Omega, X0, Y0, T, par);
     disp('Warm-up finished.');
     
     ID1 = tic;
     par.Max_iter = 200;
     [~, ~, ~, ~, Lambda1, Res1, ~, time_mat] = LOBPCG_C(Omega, X0, Y0, T, par);
     time(1, i) = toc(ID1)
     C_mat = time_mat/time(1, i)
     iter_C = length(Res1);
     res_C = Res1(end);
     Res_C(i, 1 : iter_C) = Res1;
    
    %% LOBPCG-OmegaIHL
    name = "LOBPCG-OmegaIHL";
    disp(name);
    par.reorth_Omega = 1; % the number of reorthogonalization
    par.Max_iter = 2;
    disp('Warming up...');
    LOBPCG_OmegaIHL(Omega, X0, Y0, T, par);
    disp('Warm-up finished.');

    ID2 = tic;
    par.Max_iter = 200;
    [~, ~, ~, ~, Lambda2, Res2, ~, time_mat] = LOBPCG_OmegaIHL(Omega, X0, Y0, T, par);
    time(2, i) = toc(ID2)
    OmegaIHL_mat = time_mat/time(2, i)
    iter_OmegaIHL = length(Res2);
    res_OmegaIHL = Res2(end);
    Res_OmegaIHL(i, 1 : iter_OmegaIHL) = Res2;

    %% LOBPCG-CIHL
    name = "LOBPCG-CIHL";
    disp(name);
    par.Max_iter = 2;
    disp('Warming up...');
    LOBPCG_CIHL(Omega, X0, Y0, T, par);
    disp('Warm-up finished.');
    
    ID3 = tic;
    par.Max_iter = 200;
    [~, ~, ~, ~, Lambda3, Res3, ~, time_mat] = LOBPCG_CIHL(Omega, X0, Y0, T, par);
    time(3, i) = toc(ID3)
    CIHL_mat = time_mat/time(3, i)
    iter_CIHL = length(Res3);
    res_CIHL = Res3(end);
    Res_CIHL(i, 1 : iter_CIHL) = Res3;
    
    %% LOBPCG-CIHL-OmegaIHL
    name = "LOBPCG-CIHL-OmegaIHL";
    disp(name);
    par.reorth_Omega = 1; % the number of reorthogonalization
    par.Max_iter = 2;
    par.criterion = 4; % turning algorithm criterion
    disp('Warming up...');
    LOBPCG_CIHL_OmegaIHL(Omega, X0, Y0, T, par);
    disp('Warm-up finished.');

    ID4 = tic;
    par.Max_iter = 200;
    [~, ~, ~, ~, Theta4, Res4, ~, time_mat, arg_CIHL_OmegaIHL] = LOBPCG_CIHL_OmegaIHL(Omega, X0, Y0, T, par);
    time(4, i) = toc(ID4)
    CIHL_OmegaIHL_mat = time_mat/time(4, i)
    iter_CIHL_OmegaIHL = length(Res4);
    res_CIHL_OmegaIHL = Res4(end);
    Res_CIHL_OmegaIHL(i, 1: iter_CIHL_OmegaIHL) = Res4;
    
    %% figure
    figure(i)
    semilogy(Res1, 'r-', 'Linewidth', 2, 'DisplayName', 'LOBPCG-C'); % 红色实线
    hold on;
    semilogy(Res2, 'm--', 'Linewidth', 2, 'DisplayName', 'LOBPCG-OmegaIHL'); %品红色虚线
    hold on;
    semilogy(Res3, 'b-.', 'Linewidth', 2, 'DisplayName', 'LOBPCG-CIHL'); % 蓝色点划线
    hold on;
    semilogy(Res4, '-', 'Linewidth', 2, 'DisplayName', 'LOBPCG-CIHL-OmegaIHL'); % 黑色点线

    xlabel('Iteration');
    ylabel('Residual');

    legend show;
    legend('FontSize', 8);
    
    
    iter_C
    iter_OmegaIHL
    iter_CIHL
    arg_CIHL_OmegaIHL
    
    
    nk = [n, par.l_pos];
    Iteration = [Iteration; [iter_C, iter_OmegaIHL, iter_CIHL, iter_CIHL_OmegaIHL]];
    Iter_Res_Time = [iter_C, iter_OmegaIHL, iter_CIHL, iter_CIHL_OmegaIHL;
                      res_C,  res_OmegaIHL, res_CIHL, res_CIHL_OmegaIHL;
                      time(1, i), time(2, i), time(3, i), time(4, i)];
    time_mat = [C_mat, OmegaIHL_mat, CIHL_mat, CIHL_OmegaIHL_mat];
    
    
    dlmwrite('./Data/LOBPCG_iter_res_time.txt', nk, '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite('./Data/LOBPCG_iter_res_time.txt', Iter_Res_Time, '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite('./Data/LOBPCG_iter_res_time.txt', time_mat, '-append', 'delimiter', ',', 'precision', 4);
    
    Res = sprintf('./Data/Res_case%d.txt', i);
    Iter = sprintf('./Data/Iter_case%d.txt', i);
    dlmwrite(Iter, Iteration(i, :), '-append', 'delimiter', ',');
    dlmwrite(Res, Res_C(i, :), '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite(Res, Res_OmegaIHL(i, :), '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite(Res, Res_CIHL(i, :), '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite(Res, Res_CIHL_OmegaIHL(i, :), '-append', 'delimiter', ',', 'precision', 4);
    
end

hold off;
