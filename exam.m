% given information
mass = [60.1429 18.5714;
    18.5714 14.2857];

damping = [0 0;
    0 2];

spring = [0 0;
    0 1000];

curvy_B = [1; 0];
curvy_F = [-.042857; -.071429];

% state space matrices
A = [zeros(2,2) eye(2,2);
    -inv(mass)*spring -inv(mass)*damping];

B = [zeros(2,2);
    mass\curvy_B mass\curvy_F];

M = [1 0 0 0];

base_sys = ss(A, B, M, zeros(1,2));

% for parts 4, 5 and 6, f can be ignored
B = B * [1; 0]; % remove effects of f from B

% given desired information (we're group 4)
ts_desired = 30; % seconds
os_desired = 16; % percent
torque_limit = 3; % no units
setpoint = deg2rad(60); %radians

damp_desired = sqrt((log(os_desired/100)^2)/(pi^2 + log(os_desired/100)^2));

%% 4a: Full state feedback controller
% FSF requires finding damping ratio, natural frequency of base system
% (obtained from eigenvalues of A)
w_n = 10.8141; % rad/s
damp = 0.01081;

% with given desired settling time and %OS/damping ratio, the 2nd-order
% dominant natural frequency can be found through 
sigma = -log(.02 * sqrt(1-damp_desired^2)) / ts_desired;
wn_desired = sigma / damp_desired;

% desired dominant poles: s = (damp * w_n) +/- w_n * (sqrt(damp^2 - 1))
fsf_poles = [-damp_desired * wn_desired + wn_desired * sqrt(damp_desired^2 -1) ... 
    -damp_desired * wn_desired - wn_desired * sqrt(damp_desired^2 -1)];

% add two more poles with the same frequencies but much larger magnitudes
lambda = [-1.8939 + wn_desired * sqrt(damp_desired^2 -1) ...
    -1.8939 - wn_desired * sqrt(damp_desired^2 -1) ...
    fsf_poles];

% G: the gain matrix which places the closed-loop poles at the desired
% locations
G = place(A, B, lambda);


%% 4b: estimator for above FSF control

% instead of A-GB, it's A-FM, where F is a gain matrix like G but for M
% instead of B
F = place(A.', M.', (3 * real(lambda) + imag(lambda)*1i)).';

% estimator error as an LTI system
estimator_error = ss((A - F*M), zeros(4,1), eye(4), zeros(1,1));

% arbitrary initial conditions, with closed loop (8 state) companion
inits = [deg2rad(5), deg2rad(10), deg2rad(3), deg2rad(4)];
inits_cl = [inits inits];

[y_e, t_e] = initial(estimator_error, inits);

figure()
hold on
grid on

subplot(2,2,1)
plot(t_e, rad2deg(y_e(:,1)))
title("Torso Angle $\theta$", 'Interpreter', 'latex')
xlabel("Time(s)")
ylabel("Angle (deg)")

subplot(2,2,2)
plot(t_e, rad2deg(y_e(:,2)))
title("Hip Angle $\psi$", 'Interpreter', 'latex')
xlabel("Time(s)")
ylabel("Angle (deg)")

subplot(2,2,3)
plot(t_e, rad2deg(y_e(:,3)))
title("Torso Velocity $\dot{\theta}$", 'Interpreter', 'latex')
xlabel("Time(s)")
ylabel("Angular Velocity (deg/s)")

subplot(2,2,4)
plot(t_e, rad2deg(y_e(:,4)))
title("Hip Velocity $\dot{\psi}$", 'Interpreter', 'latex')
xlabel("Time(s)")
ylabel("Angular Velocity (deg/s)")

sgtitle("Estimator Error, All States")
hold off

%% 4c: closed loop system and set-point control

% closed-loop matrices for the ss() commands
A_cl = [A -B*G;
    F*M A-B*G-F*M];

B_cl = [B;
    B];

C_cl = [M zeros(1,4)];

sys_cl = ss(A_cl, B_cl, C_cl, zeros(1,1));

% set-point controller (sorta?)
opt = stepDataOptions('StepAmplitude', .1413);


%% 4d: zero initial conditions, estimator only

% is this just totally redundant with 5?

%% 5: zero initial conditions, system
step(sys_cl, opt)

%% 6: look at all of the outputs

M_full = [eye(4) zeros(4); zeros(4,8)];

sys_full = ss(A_cl, B_cl, M_full, zeros(1,1));

[y, t] = step(sys_full, opt);

figure()
hold on

subplot(2,2,1)
plot(t, rad2deg(y(:,1)))
title("Torso Angle $\theta$", 'Interpreter', 'latex')
xlabel("Time(s)")
ylabel("Angle (deg)")
grid on

subplot(2,2,2)
plot(t, rad2deg(y(:,2)))
title("Hip Angle $\psi$", 'Interpreter', 'latex')
xlabel("Time(s)")
ylabel("Angle (deg)")
grid on

subplot(2,2,3)
plot(t, rad2deg(y(:,3)))
title("Torso Velocity $\dot{\theta}$", 'Interpreter', 'latex')
xlabel("Time(s)")
ylabel("Angular Velocity (deg/s)")
grid on

subplot(2,2,4)
plot(t, rad2deg(y(:,4)))
title("Hip Velocity $\dot{\psi}$", 'Interpreter', 'latex')
xlabel("Time(s)")
ylabel("Angular Velocity (deg/s)")
grid on

sgtitle("System Response, All States")
hold off