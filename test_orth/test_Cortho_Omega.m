function [ortho_V, ortho_V_rel] = test_Cortho_Omega(VX, VY, Omega)
%% test Omega-orthogonal
q = size(VX, 2);
V = [VX, conj(VY);
     VY, conj(VX)];
ortho_V = norm(V'*Omega*V - eye(2*q));
ortho_V_rel = ortho_V/(norm(V)^2*norm(Omega) + 1);