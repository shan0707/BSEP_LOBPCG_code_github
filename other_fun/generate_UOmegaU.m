function [UOmegaU_A, UOmegaU_B, time_mat] = generate_UOmegaU(UX, UY, OmegaZX,...
                                OmegaZY, A, B, conj_A, conj_B, k_pos, np, time_mat)
%%[UOmegaU_A, UOmgeaUB; -conj(UOmgeaUB), -conj(UOmgeaUA)] = U'OmegaU
ZX = UX(:, 1 : k_pos);
ZY = UY(:, 1 : k_pos);

WX = UX(:, k_pos + np + 1 : end);
WY = UY(:, k_pos + np + 1 : end);

id1 = tic;
OmegaWX = A*WX + B*WY; 
OmegaWY = conj_B*WX + conj_A*WY;
time_mat = time_mat + toc(id1);

if np ~= 0
    PX = UX(:, k_pos + 1 : k_pos + np);
    PY = UY(:, k_pos + 1 : k_pos + np);
    id2 = tic;
    OmegaPX = A*PX + B*PY;
    OmegaPY = conj_B*PX + conj_A*PY;
    time_mat = time_mat + toc(id2);
end

ZOmegaZ1 = ZX'*OmegaZX + ZY'*OmegaZY;
ZOmegaZ1 = (ZOmegaZ1 + ZOmegaZ1')/2;
ZOmegaW1 = ZX'*OmegaWX + ZY'*OmegaWY;
WOmegaW1 = WX'*OmegaWX + WY'*OmegaWY;
WOmegaW1 = (WOmegaW1 + WOmegaW1')/2;

ZOmegaZ2 = ZY.'*OmegaZX + ZX.'*OmegaZY;
ZOmegaZ2 = (ZOmegaZ2 + ZOmegaZ2.')/2;
ZOmegaW2 = ZY.'*OmegaWX + ZX.'*OmegaWY;
WOmegaW2 = WY.'*OmegaWX + WX.'*OmegaWY;
WOmegaW2 = (WOmegaW2 + WOmegaW2.')/2;

%WOmegaZ1 = WX'*OmegaZX + WY'*OmegaZY;
%WOmegaZ2 = WY.'*OmegaZX + WX.'*OmegaZY;

if np == 0 % U = [Z, W]
    UOmegaU_A = [ZOmegaZ1, ZOmegaW1;
                 ZOmegaW1', WOmegaW1];
    UOmegaU_B = [ZOmegaZ2, ZOmegaW2;
                 ZOmegaW2.', WOmegaW2];
    UOmegaU_B = conj(UOmegaU_B);
else  % U = [Z, W, P]
    ZOmegaP1 = ZX'*OmegaPX + ZY'*OmegaPY;
    POmegaP1 = PX'*OmegaPX + PY'*OmegaPY;
    POmegaP1 = (POmegaP1 + POmegaP1')/2;
    POmegaW1 = PX'*OmegaWX + PY'*OmegaWY;

    ZOmegaP2 = ZY.'*OmegaPX + ZX.'*OmegaPY;
    POmegaP2 = PY.'*OmegaPX + PX.'*OmegaPY;
    POmegaP2 = (POmegaP2 + POmegaP2.')/2;
    POmegaW2 = PY.'*OmegaWX + PX.'*OmegaWY;

    %WOmegaP1 = WX'*OmegaPX + WY'*OmegaPY;
    %POmegaZ1 = PX'*OmegaZX + PY'*OmegaZY;
    %WOmegaP2 = WY.'*OmegaPX + WX.'*OmegaPY;
    %POmegaZ2 = PY.'*OmegaZX + PX.'*OmegaZY;

    UOmegaU_A = [ZOmegaZ1,  ZOmegaP1,  ZOmegaW1;
                 ZOmegaP1', POmegaP1,  POmegaW1;
                 ZOmegaW1', POmegaW1', WOmegaW1];
    %UOmegaU_A = [diag(diag(ZOmegaZ1)), zeros(k_pos, k_pos), ZOmegaW1;
    %             zeros(k_pos, k_pos), POmegaP1, POmegaW1;
    %             ZOmegaW1', POmegaW1', WOmegaW1];
    UOmegaU_B = [ZOmegaZ2,   ZOmegaP2,   ZOmegaW2;
                 ZOmegaP2.', POmegaP2,   POmegaW2;
                 ZOmegaW2.', POmegaW2.', WOmegaW2];
    %UOmegaU_B = [zeros(k_pos, k_pos), zeros(k_pos, k_pos), ZOmegaW2;
    %             zeros(k_pos, k_pos), POmegaP2, POmegaW2;
    %             ZOmegaW2.', POmegaW2.', WOmegaW2];
    
    UOmegaU_B = conj(UOmegaU_B);
end
