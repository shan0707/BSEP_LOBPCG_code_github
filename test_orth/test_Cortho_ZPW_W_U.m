function [ortho_ZPW, ortho_W, ortho_U] = test_Cortho_ZPW_W_U(ZPX, ZPY, WX, WY, UX, UY, n)
%% test orthogonal Q, Z, P, PZ, PP
ZP = [ZPX, conj(ZPY); 
      ZPY, conj(ZPX)];
W = [WX, conj(WY);  
     WY, conj(WX)];
Cn = [eye(n), zeros(n); 
     zeros(n), -eye(n)];
U = [UX, conj(UY); 
     UY, conj(UX)];
p = size(W, 2);
q = size(U, 2);
S = U'*Cn*U;
ortho_ZPW = norm(ZP'*Cn*W);
ortho_W = norm(W'*Cn*W - [eye(p/2), zeros(p/2); zeros(p/2), -eye(p/2)]);
ortho_U = norm(S - [eye(q/2), zeros(q/2); zeros(q/2), -eye(q/2)]);
