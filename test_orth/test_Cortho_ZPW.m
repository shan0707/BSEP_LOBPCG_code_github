function [ortho_ZPW, ortho_ZPW_real] = test_Cortho_ZPW(ZPX, ZPY, WX, WY)
%% test orthogonal ZP, W
n = size(WX, 1);
ZP = [ZPX, conj(ZPY); 
      ZPY, conj(ZPX)];
W = [WX, conj(WY);  
     WY, conj(WX)];

Cn = [speye(n), sparse(n, n); 
     sparse(n, n), -speye(n, n)];
ortho_ZPW = norm(ZP'*Cn*W);
ortho_ZPW_real = norm(ZP'*Cn*W)/(norm(ZP)*norm(W));