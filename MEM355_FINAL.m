k=[0 0; 0 1000];
d=[0 0; 0 2];
T=[1 0];
f=[-0.042857 -0.071429];
m=[60.1429 18.5714; 18.5714 14.2857];
aPart1=-k/m;
aPart2=-d/m;
bT= T/m;
bF= f/m;
a=[0 0 1 0; 0 0 0 1; 0 36.110 0 0.07222; 0 -116.9443 0 -0.23389];
b=[0 0; 0 0; bT; bF];
M=[1 0 0 0];
%Question 2 Code
[V,D]=eig(a);
eigenvector1=a*V;
[d, ind]= sort(diag(D));
e=eig(a);
sys = ss(a,b(:,2),eye(4),[0 0 0 0]');
%Part 4a
desiredOS=16;
desiredTsettle=30;
dP1=-0.1333333+0.2285732i
dP2=-0.1333333-0.2285732i
%part 4c and what you need to change
nP1=-1.8939+.2319i
nP2=-1.893-.2319i
G=[dP1 dP2 nP1 nP2];
L=place(a,b,G);
char = [];
eigs = eig(a-b*L);
[max_real,i] = min(eigs);
dom_pole = eigs(i);
nat_freq = abs(dom_pole);
zeta = -max_real/nat_freq;
t_const  = 1/abs(max_real);
t_settle = 4*t_const;
per_OS = 100*exp(-zeta*pi/sqrt(1-zeta^2));
char = [char;per_OS zeta nat_freq t_settle t_const];
Closed_A = a - b*L;
char(1,:)
closed_sys = ss(Closed_A,b(:,1),eye(4),[0 0 0 0]');
[y,t,x]=step(closed_sys,400);
%[y2,t2,x2]=step(sys,400)
D=rad2deg(x)
%D2=rad2deg(x2)
plot(t,D(:,1)), title('step response with controller'),ylabel('torso position (degrees)'); grid
%hold on
%plot(t2,D2(:,1))
%hold off
%subplot(4,1,1), plot(t,D(:,1)), title('step thrust response'),ylabel('torso position (deg)')
%subplot(4,1,2), plot(t,D(:,2)),ylabel('leg position (deg)')
%subplot(4,1,3), plot(t,D(:,3)),ylabel('torso velocity (deg/sec)')
%subplot(4,1,4), plot(t,D(:,4)),ylabel('leg velocity (deg/sec)'),
%xlabel('time(sec)')