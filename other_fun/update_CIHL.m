function [ZX, ZY, PX, PY] = update_CIHL(UX, UY, VX, VY, k_pos)
%% Update P, Z by the improved HL trick
index1 = 1 : k_pos;
index2 = (k_pos + 1) : size(VX, 2);

VX1 = VX(:, index1);
VX2 = VX(:, index2);

VY1 = VY(:, index1);
VY2 = VY(:, index2);

conj_UX = conj(UX);
conj_UY = conj(UY);

%% update Z
ZX = UX*VX1 + conj_UY*VY1;
ZY = UY*VX1 + conj_UX*VY1;

%% update Q
VX12 = VX2(index1, :)';
VY12 = -VY2(index1, :).';

l2 = length(index2);
if l2 < length(index1)
    VX12 = VX12(:, 1 : l2);
    VY12 = VY12(:, 1 : l2);
end

[QX, QY] = ortho_C_SVQB(VX12, VY12);
[QX, QY] = ortho_C_SVQB(QX, QY);

%% update PP
PPX = VX2*QX + conj(VY2)*QY;
PPY = VY2*QX + conj(VX2)*QY;

%% update P
PX = UX*PPX + conj_UY*PPY;
PY = UY*PPX + conj_UX*PPY;