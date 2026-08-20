function [max_res, max_res_locked, locked_id] = compute_Res(Omega, Omega_norm, V, Lambda, n, l_pos, tol)
%% Omega*V-Cn*V*Lambda
CV = [V(1 : n, :); - V(n + 1 : 2*n, :)];
locked_id = [];
max_res_locked = 0;
for i = 1 : l_pos
    R = Omega*V(:, i) - CV(:, i)*Lambda(i);
    res_norm = norm(R)/((Omega_norm + abs(Lambda(i)))*norm(V(:, i)));
    if res_norm < tol
        if i == 1
            locked_id = 1;
        else
            locked_id = [locked_id; i];
        end
        max_res_locked = max(max_res_locked, res_norm);
    end
    if i == 1
        max_res = res_norm;
    end
    if res_norm > max_res
        max_res = res_norm;
    end
end