function [ortho_Q, ortho_Z, ortho_P, ortho_PZ, ortho_PP] = test_Cortho_Q_P_Z...
                                        (QX, QY, PX, PY, ZX, ZY, PPX, PPY)
%% test orthogonal Q, Z, P, PZ, PP
Q = [QX, conj(QY); 
     QY, conj(QX)];
p = size(Q, 1)/2;
q = size(Q, 2)/2;
Cp = [eye(p), zeros(p);
      zeros(p), -eye(p)];
Cq = [eye(q), zeros(q);
      zeros(q), -eye(q)];
ortho_Q = norm(Q'*Cp*Q - Cq);

% test P
P = [PX, conj(PY); 
     PY, conj(PX)];
p = size(P, 2)/2;
n = size(P, 1)/2;
Cn = [eye(n), zeros(n);
      zeros(n), -eye(n)];
Cp = [eye(p), zeros(p);
      zeros(p), -eye(p)];
ortho_P = norm(P'*Cn*P - Cp);

% test PP
PP = [PPX, conj(PPY); 
      PPY, conj(PPX)];
pp = size(PP, 2)/2;
n = size(PP, 1)/2;
Cn = [eye(n), zeros(n);
      zeros(n), -eye(n)];
Cpp = [eye(pp), zeros(pp);
      zeros(pp), -eye(pp)];
ortho_PP = norm(PP'*Cn*PP - Cpp);

%test Z
Z = [ZX, conj(ZY); 
     ZY, conj(ZX)];
z = size(Z, 2)/2;
n = size(Z, 1)/2;
Cn = [eye(n), zeros(n);
      zeros(n), -eye(n)];
Cz = [eye(z), zeros(z);
      zeros(z), -eye(z)];
ortho_Z = norm(Z'*Cn*Z - Cz);

%test P, Z
ortho_PZ = norm(P'*Cn*Z);