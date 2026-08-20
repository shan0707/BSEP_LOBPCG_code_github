function [UX, UY] = ortho_Omega_CGS(UX, UY, Omega, reorth)
%% Omega-orthonormalize [UX, conj(UY); UY, conj(UX)] use the CGS
[n, p] = size(UX);

OmegaUX = zeros(n, p);
OmegaUY = zeros(n, p);

steps = 1;
if nargin >= 4
    steps = 1 + reorth;
end
OmegaU = Omega*[UX; UY];

for i = 1 : p
    ui = [UX(:, i); UY(:, i)];
    Omegaui = OmegaU(:, i);

    if i > 1
        id = 1 : i - 1;
        UX_id = UX(:, id);
        UY_id = UY(:, id);
        OmegaUX_id = OmegaUX(:, id);
        OmegaUY_id = OmegaUY(:, id);

        U_id = [UX_id, conj(UY_id);
                UY_id, conj(UX_id)];

        OmegaUid = [OmegaUX_id, conj(OmegaUY_id);
                    OmegaUY_id, conj(OmegaUX_id)];

        for k = 1 : steps
            c = U_id'*Omegaui;
            ui = ui - U_id*c;
            Omegaui = Omegaui - OmegaUid*c;
        end
    end

    % Omega-orthonormalize vectors and maintain the required symmetry structure.
    scale = sqrt(abs(ui'*Omegaui)); 
    ui = ui/scale;
    Omegaui = Omegaui/scale;

    x = ui(1 : n);       
    y = ui(n + 1 : end);
    Omegax = Omegaui(1 : n);    
    Omegay = Omegaui(n + 1 : end);

    u = [x, conj(y);
         y, conj(x)];
    Omegau = [Omegax, conj(Omegay);
              Omegay, conj(Omegax)];

    uOmegau = u'*Omegau; 
    [z, ~] = mysqrtm(uOmegau);

    u = u/z;
    Omegau = Omegau/z;

    UX(:, i) = u(1 : n, 1);
    UY(:, i) = u(n + 1 : end, 1);
    OmegaUX(:, i) = Omegau(1 : n, 1);
    OmegaUY(:, i) = Omegau(n + 1 : end, 1);
end
end
