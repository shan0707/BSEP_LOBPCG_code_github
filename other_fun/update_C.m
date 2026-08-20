function [ZX, ZY, PX, PY] = update_C(UX, UY, VX, VY, k_pos)
%% Update P, Z
index1 = 1 : k_pos;
index2 = (k_pos + 1) : size(VX, 1);

VX1 = VX(index1, :);
VX2 = VX(index2, :);

VY1 = VY(index1, :);
VY2 = VY(index2, :);

UX1 = UX(:, index1);
UX2 = UX(:, index2);

UY1 = UY(:, index1);
UY2 = UY(:, index2);

conj_UX1 = conj(UX1);
conj_UX2 = conj(UX2);

conj_UY1 = conj(UY1);
conj_UY2 = conj(UY2);

PX = UX2*VX2 + conj_UY2*VY2;
PY = UY2*VX2 + conj_UX2*VY2;
ZX = UX1*VX1 + conj_UY1*VY1 + PX;
ZY = UY1*VX1 + conj_UX1*VY1 + PY;