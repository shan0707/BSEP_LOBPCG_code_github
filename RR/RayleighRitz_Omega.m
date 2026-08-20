function [Theta, VX, VY] = RayleighRitz_Omega(UX, UY, k_pos)
%% Rayleigh-Ritz with Omega-inner product 
%output the k_pos largest eigenpair

CX = UX'*UX - UY'*UY;
CX = (CX + CX')/2;
T = UX'*conj(UY);
CY = T - T.';

if nargin > 2
    [VX, VY, Theta] = Similar_BSE_complex(CX, CY, k_pos);  % output k_pos positive eigenpairs
else
    [VX, VY, Theta] = Similar_BSE_complex(CX, CY); % output all positive eigenpairs
end