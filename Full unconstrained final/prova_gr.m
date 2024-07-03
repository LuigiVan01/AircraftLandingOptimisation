clear all
close all
clc

%% Model parameters
M = 22e3;
m = 130;
J = 100e3;
k = 6.73e5;
k_w = 1.59e6;
c = 4066;
c_w = 1.43e5;
L = 10;
S = 1/2*10^2;
Lf = 7.76;
Lr = 1.94;
T_max = 6e4;
theta_max = 20/180*pi;
Brake_max = 2e5;
Fa_max = 1e5;

th = [M J m k k_w c c_w S L Lr Lf T_max theta_max Brake_max Fa_max]';
g = 9.81;

d = 0;

%% Control problem parameters
Ts          =   0.01;                            % Sampling time
Tend        =   15;                               % Time horizon
Tend_td     =   1;
Tend_b      =   14;
N           =   Tend/Ts;                         % Prediction steps
N_td        =   Tend_td/Ts;
N_b         =   Tend_b/Ts;
ds_u        =   50;
ds_u_brake        =  70;
ds_u_touchdown    =  50;
nz          =   6;
nu          =   6;
n_free      =   4;
x_ref       =   600;
zdd_w  =   0;
thdd_w =   0*180/pi;
Q           =   diag([0;0;0;zdd_w;0;thdd_w]);
Qf_dot_gr       =   diag([1e3;1e3;0;0;0;0]);
Qf_gr           =   diag([0;0;0;0;0;0]);       % Terminal weight
R           =   zeros(nu*N/ds_u);                  % Input weight

z = zeros(nz,N+1);
zd = zeros(nz,N+1);
z0_free          = [95; -0.5; 5/180*pi; -2/180*pi];
Xstar(1:n_free,1) = z0_free;
z0_star = [0;Xstar(1,1);0;Xstar(2,1);Xstar(3,1);Xstar(4,1)];
z(:,1) = z0_star;
i = 1;
time = 0:Ts:Tend;
Nu = nu*(N_td/ds_u_touchdown + N_b/ds_u_brake);
u = zeros(nu,N);
ztemp = z0_star;

u(1,:) = 0;
u(2,:) = 0.1;
u(3,:) = 0.8;
u(4,:) = 0;
u(5,:) = 0;
u(6,:) = 0;

U0_tdp = zeros(nu*N_td/ds_u_touchdown,1);
U0_bp = zeros(nu*N_b/ds_u_brake,1);
for i = 1:nu
    for j = 1:N_td/ds_u_touchdown
        U0_tdp(i+(j-1)*nu,1) = u(i,j);
        if i==5 || i==6
            U0_tdp(i+(j-1)*nu,1) = 1;
        end
    end
end

for i = 1:nu
    for j = 1:N_b/ds_u_brake
        U0_bp(i+(j-1)*nu,1) = u(i,j+N_td/ds_u_touchdown);
    end
end

%% choice of the brake force
b = 0.8;

%% Check of the costs of the input sequence
% F = Ground_cost_GN(Xstar,n_free,nu,d,Ts,Tend_td,Tend_b,ds_u_touchdown,ds_u_brake,Q,R,Qf_gr,Qf_dot_gr,x_ref,th);
% F_state = F(1:nz*N,1);
% F_final_state = F(nz*N+1:nz*(N+1),1)
% F_final_dot = F(nz*(N+1)+1:nz*(N+2),1)
% F_barrier = F(nz*(N+2)+1:end,1);

%% Simulation
% RK2 simulation
u_check = zeros(nu,1);
for i = 1:N_td/ds_u_touchdown
    for ind = (i-1)*ds_u_touchdown+2:i*ds_u_touchdown+1
        if z(2,ind-1) > 0
            U0_tdp((i-1)*nu+4,1) = b;
        end
        u_check = U0_tdp((i-1)*nu+1:i*nu,1);
        z_dot = ground2(0, z(:,ind-1), u_check, d, th);
        u(:,ind-1) = u_check;
        zd(:,ind) = z_dot;
        z_prime = z(:,ind-1) + Ts/2*z_dot;
        z(:,ind) = z(:,ind-1) + Ts*ground2(0, z_prime, u_check, d, th);
    end
end

for i = 1:N_b/ds_u_brake
    for ind = N_td + (i-1)*ds_u_brake+2:N_td + i*ds_u_brake+1
        if z(2,ind-1) > 2
        U0_bp((i-1)*nu+4,1) = b;
        end
        u_check = U0_bp((i-1)*nu+1:i*nu,1);
        z_dot = ground2(0, z(:,ind-1), u_check, d, th);
        u(:,ind-1) = u_check;
        zd(:,ind) = z_dot;
        z_prime = z(:,ind-1) + Ts/2*z_dot;
        z(:,ind) = z(:,ind-1) + Ts*ground2(0, z_prime, u_check, d, th);
    end
end

%% Saving the initial guess

% Xstar = [z0_free;U0_tdp;U0_bp];
% save initialguess_gr.mat Xstar

%% Simultion data from the optimization function (to check the simulation of the script)
[~,z_sim] = Ground_cost_GN(Xstar,n_free,nu,d,Ts,Tend_td,Tend_b,ds_u_touchdown,ds_u_brake,Q,R,Qf_gr,Qf_dot_gr,x_ref,th);
zsim = zeros(nz,N+1);
for ind = 1:N+1
    zsim(:,ind) = z_sim((ind-1)*nz+1:ind*nz,1);
end

%% Plots of th simulation
for j=1:nz
    figure(j)
    if j == 5 || j == 6
        hold on
        plot(time,180/pi*z(j,:),'b');
        plot(time,180/pi*zsim(j,:),'g--o');
        plot(time,180/pi*zd(j,:),'r');
        hold off
    else
        hold on
        plot(time,z(j,:),'b');
        plot(time,zsim(j,:),'g--o');
        plot(time,zd(j,:),'r');
        hold off
    end
end
for j=1:nu
    figure(nz+j)
    plot(u(j,:));
end
figure(13)
hold on
plot(z(1,:), z(3,:),'g');
plot(zsim(1,:), zsim(3,:),'b--o');
hold off
    

