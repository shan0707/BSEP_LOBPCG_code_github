function [WX, WY] = ortho_together_C(WX, WY, ZPX, ZPY)
%% Cn-orthogonalize [WX, conj(WY); WY, conj(WX)] about [ZPX, conj(ZPY); ZPY, conj(ZPX)]
%% and Cn-orthonormalize [WX, conj(WY); WY, conj(WX)]
n = size(WX, 1);

b = randn(size(ZPX, 2), 5);

conj_ZPY = conj(ZPY);
conj_ZPX = conj(ZPX);
CnZP = [ZPX,  conj_ZPY; 
       -ZPY, -conj_ZPX];
ZPCk = [ZPX, -conj_ZPY; 
        ZPY, -conj_ZPX];

ZPXb = ZPX*b; ZPYb = ZPY*b;
ZP_norm = norm([ZPXb; ZPYb], 'fro')/norm(b, 'fro');

for i = 1 : 2
    W = [WX; WY];
    W = W - ZPCk*(CnZP'*W);
    WX = W(1 : n, :);
    WY = W(n + 1 : end, :);
   
    for j = 1 : 2
        [WX, WY] = ortho_C_SVQB(WX, WY);
        %WXa = WX*a;
        %WYa = WY*a;
        %ErrX = WX'*WXa - WY'*WYa - a;
        %ErrY = WY.'*WXa - WX.'*WYa;
        %W_norm = norm([WXa; WYa], 'fro')/norm(a, 'fro');
        %ortho_W = (norm(ErrX, 'inf') + norm(ErrY, 'inf'))/(W_norm^2);
        %if ortho_W < 1e-15
        %    break;
        %end
    end
    
    ErrX = WX'*ZPXb - WY'*ZPYb;
    ErrY = WY.'*ZPXb - WX.'*ZPYb;
    probe = randn(size(WX, 2), 5);
    W_norm = norm([WX; WY]*probe, 'fro')/norm(probe, 'fro');
    ortho_WZP = (norm(ErrX, 'inf') + norm(ErrY, 'inf'))/(W_norm*ZP_norm);

    if ortho_WZP < 1e-15
         break;
    end
end
