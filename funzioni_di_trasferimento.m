clear all
close all
clc

i=sqrt(-1);
s=tf('s');

% G11=(-0.045)/(s^2+0.33*s+1.5);
% G12=(0.024)/(s^2+0.33*s+1.5);
% G21=0.04/(14*s+1);
% G22=0.20/(14*s+1);
%  G=[G11 G12;
%      G21/s G22/s];
% sys=ss(G);

n=4;
p=2;
m=2;
A=[0 0 1 0;
   0 0 0 1;
   -1.5 0 -0.33 0;
   0 0 0 -1/14];
B=[0 0;
    0 0;
    -0.045 0.024;
    0.04/14 0.20/14];
C=[1 0 0 0;
    0 1 0 0];
D=[0 0;
    0 0];
sys=ss(A,B,C,D);

% %%PP_NO_obs
% n=4;
% m=2;
% p=2;
% Mreach=ctrb(A,B);
% rank(Mreach)==n
% %ans logical = 1 means reachable
% Mobs=obsv(A,C);
% rank(Mobs)==n
% %ans logical = 1 means observable
% pp_desired_poles=[-1.4, -1.4, -1.44, -1.44];
% Kp=place(A, B, pp_desired_poles);

% %%PP_NO_ext
% n=4;
% m=2;
% p=2;
% Mreach=ctrb(A,B);
% rank(Mreach)==n
% %ans logical = 1 means reachable
% Mobs=obsv(A,C);
% rank(Mobs)==n
% %ans logical = 1 means observable
% obsv_poles=[-80, -80, -100, -100];
% L=place(A', C', obsv_poles)';
% A_ob=A-L*C;
% B_ob=[B-L*D, L];
% C_ob=eye(n);
% D_ob=zeros(n, m+p);
% pp_desired_poles=[-1.4, -1.4, -1.44, -1.44];
% Kp=place(A, B, pp_desired_poles);

%%PP_enlarged
n=4;
m=2;
p=2;
Mreach=ctrb(A,B);
rank(Mreach)==n
%ans logical = 1 means reachable
Mobs=obsv(A,C);
rank(Mobs)==n
%ans logical = 1 means observable
% obsv_poles=[-83.362796000979470 + 83.362848950906010i, -83.362796000979470 - 83.362848950906010i, -1.599768573804904e+02 + 1.599813397500049e+02i, -1.599768573804904e+02 - 1.599813397500049e+02i]; 
obsv_poles=[-80, -80, -100, -100]; %no_SAT
L=place(A', C', obsv_poles)';
A_ob=A-L*C;
B_ob=[B-L*D, L];
C_ob=eye(n);
D_ob=zeros(n, m+p);
% pp_desired_poles=[-3.2757 + 1.2122i, -3.2757 - 1.2122i, -3.3656 + 0.0078i, -3.3656 - 0.0078i]; %SAT 
pp_desired_poles=[-1.4, -1.4, -1.44, -1.44]; %no_SAT
Kp=place(A, B, pp_desired_poles);
%Enlarged system with two integrators
A_tilde=[A, zeros(n,p);
        -C, zeros(p,p)];
B_tilde=[B;
    -D];
M_tilde=[zeros(n,p);
    eye(p)];
en_desired_poles=[pp_desired_poles, -1.5, -1.5]; %[-1.5, -1.5] no_SAT; [-3.4, -3.4044] SAT
% %now we try with one integraton only, on pitch velocity [FAILURE]
% A_tilde=[A, zeros(n,1);
%     -C(1,:), zeros(1,1)];
% B_tilde=[B;
%     0 0];
% M_tilde=[zeros(n,1);
%     1; 0];
% en_desired_poles=[pp_desired_poles, -2.002];
Ken=place(A_tilde, B_tilde, en_desired_poles);
Ken_x=Ken(:, 1:n);
Ken_eta=Ken(:, n+1:end);

%%PP_with_RED_ORD_obs
T=eye(n);
At=A;
Bt=B;
Ct=C;
At11=At(1:p, 1:p);
At12=At(1:p, p+1:end);
At21=At(p+1:end, 1:p);
At22=At(p+1:end, p+1:end);
Bt1=Bt(1:p, :);
Bt2=Bt(p+1:end, :);
red_desired_poles=[-40, -40]; %[-40, -40] no_SAT; [-50, -50] SAT
Lt=place(At22', At12', red_desired_poles)';
A_xi=At22-Lt*At12;
B_xi=[Bt2-Lt*Bt1, At21+(At22-Lt*At12)*Lt-Lt*At11];
sys_red=ss(A_xi, B_xi, eye(2), zeros(2,4)); %to check if correct orderD_ob

% %%LQG_with_filter
% rank(ctrb(A, B)) == n
% % Q = eye(n); %n+p due to the enlarged system
% Q = diag([80 50 100 100]);
% % R = eye(m);
% R = diag([1 1]);
% rank(obsv(A, sqrt(Q))) == n
% Klq = lqr(A, B, Q, R);
% Klq_x = Klq(:, 1:n);
% rank(obsv(A, C)) == n
% eig(A-B*Klq)

% %% LQG_NO_LTR
% rank(ctrb(A_tilde + 0*eye(n+p), B_tilde)) == n+p 
% Q = eye(n+p);
% R = eye(m);
% rank(obsv(A_tilde + 0*eye(n+p), sqrt(Q))) == n+p
% Klq = lqr(A_tilde + 0*eye(n+p), B_tilde, Q, R);
% Klq_x = Klq(:, 1:n);
% Klq_eta = Klq(:, n+1:end);
% rank(obsv(A, C)) == n
% Q_tilde = eye(n);
% R_tilde = eye(p);
% rank(ctrb(A, sqrt(Q_tilde))) == n
% Lkf = lqr(A.', C.', Q_tilde, R_tilde).';
% A_kf = A - Lkf*C;
% B_kf = [ B - Lkf*D, Lkf];
% C_kf = eye(n);
% D_kf = zeros(n, m+p);

%% LQG_with_LTR (1.7*eye SAT; 0.73*eye no_SAT)
rank(ctrb(A_tilde + 0.73*eye(n+p), B_tilde)) == n+p 
%Q = eye(n+p); %n+p due to the enlarged system
% Q = diag([80, 50, 100, 100, 1, 1]); %SAT
Q=diag([200, 40, 1, 1, 1, 1]); %no_SAT
%R = eye(m);
R = diag([1 1]);
rank(obsv(A_tilde + 0.73*eye(n+p), sqrt(Q))) == n+p
Klq = lqr(A_tilde + 0.73*eye(n+p), B_tilde, Q, R);
Klq_x = Klq(:, 1:n);
Klq_eta = Klq(:, n+1:end);
rank(obsv(A, C)) == n
% Q_tilde = eye(n); %NO LTR
alpha=1000000000000; %1000000000000 no_SAT; 2000000000000 SAT
Q_tilde_LTR=alpha*B*B';
R_tilde = 14*eye(p); %14 no_SAT; 2 SAT
% rank(ctrb(A, sqrt(Q_tilde))) == n %NO LTR
rank(ctrb(A, sqrt(Q_tilde_LTR))) == n
% Lkf = lqr(A.', C.', Q_tilde, R_tilde).'; %NO LTR
% rho=100;
% Lkf_LTR = rho*B;
% A_kf = A - Lkf_LTR*C;
% B_kf = [ B - Lkf_LTR*D, Lkf_LTR];
Lkf = lqr(A.', C.', Q_tilde_LTR, R_tilde).';
A_kf = A-Lkf*C;
B_kf = [B-Lkf*D, Lkf];
C_kf = eye(n);
D_kf = zeros(n, m+p); %Questo è un commento
poles_sys=eig(A_tilde-B_tilde*Klq);
poles_kf=eig(A-Lkf*C);
