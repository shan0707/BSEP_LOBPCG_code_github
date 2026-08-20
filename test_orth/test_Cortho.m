function [ortho_V, ortho_V_rel] = test_Cortho(VX, VY)
%% test Cn-orthogonal
q = size(VX, 2);
ortho_V_rel = 1;
%% way 1 exact
% [p, q] = size(VX);
% Cp = [speye(p), zeros(p); 
%      zeros(p), -speye(p)];
% Cq = [speye(q), zeros(q); 
%       zeros(q), -speye(q)];
% V = [VX, conj(VY);
%      VY, conj(VX)];
% ortho_V = norm(V'*Cp*V - Cq);
% ortho_V_rel = ortho_V/(norm(V)^2);

%%
% tic
% Rand = randn(q, 7);
% VXR = VX*Rand;
% VYR = VY*Rand;
% C_I = VX'*VXR - VY'*VYR - Rand;
% C_zeros = VY.'*VXR - VX.'*VYR;
% norm_F = norm(Rand, 'fro');
% ortho_V = (norm(C_I, 'fro') + norm(C_zeros, 'fro'))/norm_F
% toc

%%
p = 1;
Rand = randn(q, p);
VXR = VX*Rand;
VYR = VY*Rand;
C_I = VX'*VXR - VY'*VYR - Rand;
C_zeros = VY.'*VXR - VX.'*VYR;
ortho_V = norm(C_I, 'inf') + norm(C_zeros, 'inf');
%V = [VX, conj(VY);
%     VY, conj(VX)];
%ortho_V_rel = ortho_V/(norm(V)^2);
