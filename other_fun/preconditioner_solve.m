function TR = preconditioner_solve(L, U, R, label)
%% solve T^{-1}R
if label == 2 && label == 4
    TR = U\(L\R);

elseif label == 3
    %incomplete cholesky decomposition
    %only for T = blkdiag(A, conj(A)); 
    n = size(R, 1)/2;
    R1 = R(1 : n, :);
    R2 = R(n + 1 : 2*n, :);
    TR1 = U\(L\R1);
    L_c = conj(L);
    U_c = conj(U);
    TR2 = U_c\(L_c\R2);
    TR = [TR1; TR2];
    
end