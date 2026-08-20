function [M, H] = generate_ex2(n)
G = zeros(n);
K = zeros(n);
INVM = diag(2*ones(1, n));
v = 0.0306;
for i = 1:n
    for j = 1:n
        if mod(i + j, 2) == 1
           G(i, j) = 4*i*j*v/(i^2-j^2);
        end
    end
end
for i = 1:n
    K(i, i) = i^2*pi^2*(1-v^2)/2;
end
H = [-1/2*G*INVM, 1/4*G*INVM*G-K;
      INVM, -1/2*INVM*G];
J = [zeros(n), eye(n);
    -eye(n), zeros(n)];
M = J*H;
M = (M + M')/2;




