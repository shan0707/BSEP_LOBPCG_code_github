function [L, U] = preconditioner_LU(T, label)
%% solve T^{-1}R
if label == 2
    % ilu
    options.type = 'ilutp';
    options.droptol = 1e-6; 
    options.udiag = true;
    options.thresh = 0;
    [L, U] = ilu(T, options);

elseif label == 3
    %incomplete cholesky decomposition
    %only for T = blkdiag(A, conj(A)); 
    n = size(T, 1)/2;
    T1 = T(1 : n, 1 : n);
    comps = [1e-3, 0.01, 0.1];
    success = false;
    for alpha = comps
         try
             opts = struct('type', 'ict', 'droptol', 1e-6, 'diagcomp', alpha);
             L = ichol(T1, opts);
             U = L';
             success = true;
             break;
         catch
             continue;
         end
     end
     if ~success
         options.type = 'ilutp';
         options.droptol = 1e-6;
         options.udiag = true;
         options.thresh = 0;
         [L, U] = ilu(T1, options);
     end
     
elseif label == 4
    %incomplete cholesky decomposition
    comps = [1e-3, 0.01, 0.1];
    success = false;
    for alpha = comps
         try
             opts = struct('type', 'ict', 'droptol', 1e-6, 'diagcomp', alpha);
             L = ichol(T, opts);
             U = L';
             success = true;
             break;
         catch
             continue;
         end
     end
     if ~success
         options.type = 'ilutp';
         options.droptol = 1e-6;
         options.udiag = true;
         options.thresh = 0;
         [L, U] = ilu(T, options);
     end
end