function [WX, WY, time_mat] = ortho_together_Omega(WX, WY, ZPX, ZPY, OmegaZPX, OmegaZPY, Omega, time_mat)
% Omega-orthogonalize [WX, conj(WY); WY, conj(WX)] about 
%[ZPX, conj(ZPY); ZPY, conj(ZPX)] and Omega-orthonormalize [WX, conj(WY); WY, conj(WX)]
n = size(WX, 1);

b = randn(size(ZPX, 2), 5);

conj_ZPY = conj(ZPY);
conj_ZPX = conj(ZPX);
ZP = [ZPX, conj_ZPY; 
      ZPY, conj_ZPX];
OmegaZP = [OmegaZPX, conj(OmegaZPY);
           OmegaZPY, conj(OmegaZPX)];
OmegaZPXb = OmegaZPX*b; OmegaZPYb = OmegaZPY*b;

OmegaZP_norm = norm([OmegaZPXb; OmegaZPYb], 'fro')/norm(b, 'fro');

for i = 1 : 2
    W = [WX; WY];
    W = W - ZP*(OmegaZP'*W);
    WX = W(1 : n, :);
    WY = W(n + 1 : end, :);

    for j = 1 : 2
        [WX, WY, time_W] = ortho_Omega_BSE(WX, WY, Omega, n);
        time_mat = time_mat + time_W;
        %WXa = WX*a;
        %WYa = WY*a;
        %AWXa = A*WXa;
        %BWYa = B*WYa;
        %conj_AWYa = conj(A)*WYa;
        %conj_BWXa = conj(B)*WXa;
        %ErrX = WX'*AWXa + WY'*conj_BWXa + WX'*BWYa+ WY'*conj_AWYa - a;
        %ErrY = WY.'*AWXa + WX.'*conj_BWXa + WY.'*BWYa + WX.'*conj_AWYa;
        %W_norm = norm([WXa; WYa], 'fro')/norm(a, 'fro');
        %ortho_W = (norm(ErrX, 'inf') + norm(ErrY, 'inf'))/(W_norm^2);
        %if ortho_W < 1e-15
        %    break;
        %end
    end

    ErrX = WX'*OmegaZPXb + WY'*OmegaZPYb; 
    ErrY = WY.'*OmegaZPXb + WX.'*OmegaZPYb;
    a = randn(size(WX, 2), 5);
    W_norm = norm([WX; WY]*a, 'fro')/norm(a, 'fro');
    ortho_WZP = (norm(ErrX, 'inf') + norm(ErrY, 'inf'))/(W_norm*OmegaZP_norm);
 
    if ortho_WZP < 1e-15
         break;
    end
end
