function [A, B, n] = generate_BSE(example)
switch example
    case 1
        A_re = h5read('naphthalene_32_C.h5', '/A_re');
        A_im = h5read('naphthalene_32_C.h5', '/A_im');
        B_re = h5read('naphthalene_32_C.h5', '/B_re');
        B_im = h5read('naphthalene_32_C.h5', '/B_im');
    case 2
        A_re = h5read('GaAs_128_C.h5', '/A_re');
        A_im = h5read('GaAs_128_C.h5', '/A_im');
        B_re = h5read('GaAs_128_C.h5', '/B_re');
        B_im = h5read('GaAs_128_C.h5', '/B_im');
    case 3
        A_re = h5read('BN1_2304_C.h5', '/A_re');
        A_im = h5read('BN1_2304_C.h5', '/A_im');
        B_re = h5read('BN1_2304_C.h5', '/B_re');
        B_im = h5read('BN1_2304_C.h5', '/B_im');
    case 4
        A_re = h5read('BN2_2304_C.h5', '/A_re');
        A_im = h5read('BN2_2304_C.h5', '/A_im');
        B_re = h5read('BN2_2304_C.h5', '/B_re');
        B_im = h5read('BN2_2304_C.h5', '/B_im'); 
    case 5
        A_3dim = h5read('PNR_10000_C.h5', '/hbse_a');
        B_3dim = h5read('PNR_10000_C.h5', '/hbse_b');
        A = squeeze(A_3dim(1, :, :)) + sqrt(-1)*squeeze(A_3dim(2, :, :));
        B = squeeze(B_3dim(1, :, :)) + sqrt(-1)*squeeze(B_3dim(2, :, :));
        n = size(A, 1);
end
if example < 5
    n =  size(A_re, 1);
    A = A_re + i * A_im;
    B = B_re + i * B_im;      
end
A = tril(A, -1) + tril(A, -1)' + real(diag(diag(A)));
B = tril(B) + tril(B, -1).';