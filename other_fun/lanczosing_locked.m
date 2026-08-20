function [a, b, V, p] = lanczosing_locked(Schol, p, Vtil, V_locked, tolr, Max_time, ak, bk)
% This is the Lanczos function producing matrix S^{-1}-orthogoanl V and the 
% element of tridiagonal elements of T such that
% -H^2V = VT + beta_{k+1}v_{k+1}e_k^T
% where H = JS^{-1}, T = diag(a) + diag(b, 1) + diag(b, -1);
% INPUT
%   Schol: the lower Cholesky factor of spd matrix S such that H = JS
%   p: the number of Lanczos steps performed 
%   Vtil: the initial guess, which can be a vector if the process starts at
%   nothing and a matrix whose last column is the initial vector and the
%   left is the already part of matrix V computed. 
%   ak, bk: the entries of matrix T which are already computed
%   length(ak) + 1 = length(bk) + 1 = size(Vtil,2);
% OUTPUT
%   a, b: the diagonal and subdiagonal of matrix T
%   V: the Lanczos vectors.
id10 = tic;
if nargin < 6 
    error('function needs at least 6 inputs')
end
if nargin < 7; ak = []; bk = []; end
if size(Vtil, 2) ~= length(ak) + 1 || size(Vtil,2) ~= length(bk) + 1
    error('check the compatibility of Vtil, ak, bk')
end
k = length(ak); % has already k step
b = [bk zeros(1, p + 1)];
a = [ak zeros(1, p)];
V = [Vtil(:, 1 : end - 1) zeros(size(Schol, 1), p + 1)];
W = [-Jmul(Schol'\(Schol\Vtil(:, 1 : end - 1))) zeros(size(Schol, 1), p)];

%% compute the norm of W_locked
t = size(V_locked, 2);
if ~isempty(V_locked)
    W_locked = -Jmul(Schol'\(Schol\V_locked));
    for i = 1 : t
        w_nrm(i) = W_locked(:, i)'*(Schol'\(Schol\W_locked(:, i)));
    end  
else
    W_locked = [];
end

%% initial
vtil = Vtil(:, end);
vtil = external_ortho(vtil);
b(k + 1) = sqrt(vtil'*(Schol'\(Schol\vtil)));
if b(k + 1) < tolr
    vtil = randn(size(vtil));
    vtil = external_ortho(vtil);
    b(k + 1) = sqrt(vtil'*(Schol'\(Schol\vtil)));
end
V(:, k + 1) = vtil/b(k + 1);

for j = 1 : p
    W(:, k + j) = -Jmul(Schol'\(Schol\V(:, k + j)));
    Sw = Schol'\(Schol\W(:, k + j));
    a(k + j) = W(:, k + j)'*Sw;
    if j == 1 && k == 0
        vtil = Jmul(Sw) - a(k + j)*V(:, k + j);
    else
        vtil = Jmul(Sw) - a(k + j)*V(:, k + j) - b(k + j)*V(:, k + j - 1);
        vtil = Sorthog(vtil, Schol, V(:, 1 : k + j), W(:, 1 : k + j), a(1 : k + j));
    end
    vtil = external_ortho(vtil);
    b(k + j + 1) = sqrt(vtil'*(Schol'\(Schol\vtil)));

    %% Breakdown check 
    % plan A
    while b(k + j + 1) < tolr
        vtil = randn(size(vtil));
        vtil = external_ortho(vtil);
        vtil = Sorthog(vtil, Schol, V(:, 1 : k + j), W(:, 1 : k + j), a(1 : k + j));
        b(k+j+1) = sqrt(vtil'*(Schol'\(Schol\vtil)));
    end
    % plan B
    %if b(k + j + 1) < tolr
    %     p = j;
    %     break;
    %end

    V(:, k + j + 1) = vtil/b(k + j + 1);
    if toc(id10) > Max_time
        break;
    end
end

%% external orthogonalization
function v_out = external_ortho(v_in)
    % v_in orthonormalzie V_locked, W_locked
    v_out = v_in;
    if isempty(V_locked)
        return; 
    end
    for i = 1 : size(V_locked, 2)
        r = W_locked(:, i)'*(Schol'\(Schol\v_out))/w_nrm(i);
        v_out = v_out - r*W_locked(:, i);

        r = V_locked(:, i)'*(Schol'\(Schol\v_out));
        v_out = v_out - r*V_locked(:, i);

        r_J = V_locked(:, i)'*Jmul(v_out); 
        v_out = v_out - r_J*W_locked(:, i);

        r_Jw = W_locked(:, i)'*Jmul(v_out);
        v_out = v_out + r_Jw*V_locked(:, i);
    end
end
end