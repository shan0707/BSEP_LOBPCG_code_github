function [max_res, lockID_pos, lockID_neg, R_norm] = deflation(R, Theta_pos, Z, ...
                                            Omega_norm, k_pos, l_pos, tol)
lockID_pos = [];
R_norm = ones(2*k_pos, 1);
R_norm(1) = norm(R(:, 1))/((Omega_norm + abs(Theta_pos(1, 1)))*norm(Z(:, 1)));
if R_norm(1) <= 0.1*tol
    lockID_pos = [lockID_pos, 1];
end
for j = 2 : k_pos
   R_norm(j) = norm(R(:, j))/((Omega_norm + abs(Theta_pos(j, j)))*norm(Z(:, j)));
    if ismember(j - 1, lockID_pos) && R_norm(j) <= 0.1*tol
        lockID_pos = [lockID_pos, j];
    end
end
lockID_neg = lockID_pos + k_pos;
R_norm(k_pos + 1 : end) = R_norm(1 : k_pos);
max_res = max(R_norm(1 : l_pos));