% jack's crap

mass = [60.1429 18.5714;
    18.5714 14.2857];

damping = [0 0;
    0 2];

spring = [0 0;
    0 1000];

curvy_B = [1; 0];
curvy_F = [-.042857; -.071429];

A = [zeros(2,2) eye(2,2);
    -inv(mass)*spring -inv(mass)*damping];

B = [zeros(2,2);
    -inv(mass)*curvy_B -inv(mass)*curvy_F];

M = [1 0 0 0];

base_sys = ss(A, B, M, zeros(1,2));

sys_noF = ss(A, B*[1; 0], M, zeros(1,1));

syms k1 k2 k3 k4 s

K = [k1 k2 k3 k4];

A_cl = A - (B*[1;0]) * K;

det(s .* eye(4,4) - A_cl);

lamb = [0 2 -.1333+.2286i -.1333-.2286i];