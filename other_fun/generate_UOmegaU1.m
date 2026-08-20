function [UOmegaU_A, UOmegaU_B] = generate_UOmegaU1(UX, UY, Theta, A, B, conj_A, conj_B, k_pos, np)
%%[UOmegaU_A, UOmgeaUB; -conj(UOmgeaUB), -conj(UOmgeaUA)] = U'OmegaU
ZX = UX(:, 1 : k_pos);
ZY = UY(:, 1 : k_pos);

WX = UX(:, k_pos + np + 1 : end);
WY = UY(:, k_pos + np + 1 : end);
OmegaWX = A*WX + B*WY; 
OmegaWY = conj_B*WX + conj_A*WY;

if np ~= 0
    PX = UX(:, k_pos + 1 : k_pos + np);
    PY = UY(:, k_pos + 1 : k_pos + np);
    OmegaPX = A*PX + B*PY;
    OmegaPY = conj_B*PX + conj_A*PY;
end

ZOmegaW1 = ZX'*OmegaWX + ZY'*OmegaWY;
WOmegaW1 = WX'*OmegaWX + WY'*OmegaWY;
WOmegaW1 = (WOmegaW1 + WOmegaW1')/2;

ZOmegaW2 = ZY.'*OmegaWX + ZX.'*OmegaWY;
WOmegaW2 = WY.'*OmegaWX + WX.'*OmegaWY;
WOmegaW2 = (WOmegaW2 + WOmegaW2.')/2;

%WOmegaZ1 = WX'*OmegaZX + WY'*OmegaZY;
%WOmegaZ2 = WY.'*OmegaZX + WX.'*OmegaZY;

if np == 0 % U = [Z, W]
    UOmegaU_A = [Theta, ZOmegaW1;
                 ZOmegaW1', WOmegaW1];
    UOmegaU_B = [zeros(k_pos, k_pos), ZOmegaW2;
                 ZOmegaW2.', WOmegaW2];
    UOmegaU_B = conj(UOmegaU_B);
else  % U = [Z, W, P]
    POmegaP1 = PX'*OmegaPX + PY'*OmegaPY;
    POmegaP1 = (POmegaP1 + POmegaP1')/2;
    POmegaW1 = PX'*OmegaWX + PY'*OmegaWY;

    POmegaP2 = PY.'*OmegaPX + PX.'*OmegaPY;
    POmegaP2 = (POmegaP2 + POmegaP2.')/2;
    POmegaW2 = PY.'*OmegaWX + PX.'*OmegaWY;

    %WOmegaP1 = WX'*OmegaPX + WY'*OmegaPY;
    %POmegaZ1 = PX'*OmegaZX + PY'*OmegaZY;
    %WOmegaP2 = WY.'*OmegaPX + WX.'*OmegaPY;
    %POmegaZ2 = PY.'*OmegaZX + PX.'*OmegaZY;

    UOmegaU_A = [Theta, zeros(k_pos, np), ZOmegaW1;
                 zeros(np, k_pos), POmegaP1, POmegaW1;
                ZOmegaW1', POmegaW1', WOmegaW1];
   
    UOmegaU_B = [zeros(k_pos, k_pos), zeros(k_pos, np), ZOmegaW2;
                 zeros(np, k_pos), POmegaP2, POmegaW2;
                 ZOmegaW2.', POmegaW2.', WOmegaW2];
    
    UOmegaU_B = conj(UOmegaU_B);
end