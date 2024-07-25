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

%% Optimization problem parameters
Ts          =   0.01;                 % Sampling Time                        
Tend_fl     =   12;                   % Flight Time                
N_fl        =   Tend_fl/Ts;           % Samples In-Flight                      
nu_fl       =   4;                    % In_flight inputs

Tend_gr     =   15;                   % Ground overall Time               
Tend_td     =   1;                    % Touchdown dynamics Time
Tend_b      =   14;                   % Brake dynamics Time
N_gr        =   Tend_gr/Ts;           % Ground samples              
N_td        =   Tend_td/Ts;           % Touchdown samples
N_b         =   Tend_b/Ts;            % Brake samples 
ds_u_b      =   70;                   % Downsampling brake inputs
ds_u_td     =   50;                   % Downsampling touchdown inputs
nu_gr       =   6;                    % Ground inputs

n_free      =   4;                    % Free variables in the Ground initial state
nz          =   6;                    % States
ds_u        =   80;                   % Overall (Flight) downsampling
Tend = Tend_fl + Tend_gr;

%% Flight initial condition definition
x_in = 0;          % Starting x_position;
xd_in = 95;        % Starting horizontal speed;
z_in = 50;         % Starting height;
zd_in = -5;        % Starting vertical speed;
th_in = 0;         % Starting pitch;
thd_in = 0;        % Starting pitch velocity;
z0_fl   =   [x_in; xd_in; z_in; zd_in; th_in; thd_in];      % Initial state In-Fight

%% States stage cost parameters
zdd_w  =   1;                                      % Vertical acceleration weight      
thdd_w =   1*180/pi;                               % Pitch acceleration weight
xdd_w  =   1;                                      % Longitudinal acceleration weight
Q_fl     =   diag([0;xdd_w;0;zdd_w;0;thdd_w]);     % In-Flight weighting Matrix   
Q_gr     =   diag([0;xdd_w;0;zdd_w;0;thdd_w]);     % Ground weighting Matrix

%% Inputs stage cost parameters
% X_fl = [U_fl];
% X_gr = [z0_gr; U_gr];

R_fl           =   zeros(nu_fl*N_fl/ds_u);                  % Input flight weight matrix (All set to zero)
for i = 1:length(R_fl)
    if mod(i,nu_fl) == 1
        R_fl(i,i) = 0;
    elseif mod(i,nu_fl) == 2
        R_fl(i,i) = 0;
    elseif mod(i,nu_fl) == 3
        R_fl(i,i) = 0;
    elseif mod(i,nu_fl) == 0
        R_fl(i,i) = 0;
    end
end

R_gr_td           =   zeros(nu_gr*N_td/ds_u_td);                  % Input touchdown ground weight matrix (All set to zero)
for i = 1:length(R_gr_td)
    if mod(i,nu_gr) == 1
        R_gr_td(i,i) = 0;
    elseif mod(i,nu_gr) == 2
        R_gr_td(i,i) = 0;
    elseif mod(i,nu_gr) == 3
        R_gr_td(i,i) = 0;
    elseif mod(i,nu_gr) == 4
        R_gr_td(i,i) = 0;
    elseif mod(i,nu_gr) == 5 
        R_gr_td(i,i) = 0;
    elseif mod(i,nu_gr) == 0
        R_gr_td(i,i) = 0;
    end
end

R_gr_b           =   zeros(nu_gr*N_b/ds_u_b);                  % Input braking ground weight matrix

for i = 1:length(R_gr_b)
    if mod(i,nu_gr) == 1
        R_gr_b(i,i) = 0;
    elseif mod(i,nu_gr) == 2
        R_gr_b(i,i) = 0;
    elseif mod(i,nu_gr) == 3
        R_gr_b(i,i) = 0;
    elseif mod(i,nu_gr) == 4
        R_gr_b(i,i) = 0;
    elseif mod(i,nu_gr) == 5 
        R_gr_b(i,i) = 0;
    elseif mod(i,nu_gr) == 0
        R_gr_b(i,i) = 0;
    end
end

%% Terminal weights
Qf_fl          =   diag([0;1e3;1e3;1e3;1e3;1e3]);     % Terminal weighting Matrix for In-Flight states        
Qf_dot_gr      =   diag([1e3;1e3;1e3;1e3;1e3;1e3]);   % Terminal weighting Matrix for final equilibrium
Qf_gr          =   diag([0;0;0;0;0;0]);               % Terminal weighting Matrix for Ground states        

%% Landing strip length
x_ref = 500;

%% Ground Optimization parameters
myoptions   =   myoptimset;
myoptions.ls_beta       = 0.3;        
% myoptions.ls_c          = 0.1;
% myoptions.tolgrad    	=	1e-6;
% myoptions.tolfun    	=	1e-18;
% myoptions.tolX       	=	1e-18;
% normal
myoptions.gradmethod    = 'CD';
myoptions.graddx        = 2^-17;
% user provided
% myoptions.gradmethod    = 'UP';
% myoptions.graddx        = eps^(1/3);
myoptions.nitermax      = 50;
myoptions.Hessmethod    = 'GN';
myoptions.BFGS_gamma  	= 1e-1;
myoptions.GN_sigma      = 1e-3;
myoptions.GN_funF       = @(X_gr)Ground_cost_GN(X_gr,n_free,nu_gr,d,Ts,Tend_td,Tend_b,ds_u_td,ds_u_b,Q_gr,R_gr_td,Qf_gr,Qf_dot_gr,x_ref,th);

%% Initial guess for the ground optimization

% Initialization of the initial guess vectors
z0_free          = [95; -0.5; 5/180*pi; -2/180*pi];     
U0_td          = ones(nu_gr*N_td/ds_u_td,1);
U0_b          = ones(nu_gr*N_b/ds_u_b,1);

% Improving initial guess
load('initialguess_gr.mat');
X0 = Xstar;

%% Generating code (mex files) for faster computation (to be ran only once whenever Tend changes)
% codegen Pendulum_cost -args {U,z0,d,Ts,Tend,Q,R,Qf,th}
% codegen Flight_cost_GN -args {U0,z0,Ts,Q,R,Qf,parameters}

%% Running the Ground optimization routine
tic
[X_grstar,fxstar_gr] = myfminunc(@(X_gr)Ground_cost(X_gr,n_free,nu_gr,d,Ts,Tend_td,Tend_b,ds_u_td,ds_u_b,Q_gr,R_gr_td,Qf_gr,Qf_dot_gr,x_ref,th),X0,myoptions);
% [X_grstar,fxstar_gr] = fminunc(@(X_gr)Ground_cost(X_gr,n_free,nu_gr,d,Ts,Tend_td,Tend_b,ds_u_touchdown,ds_u_brake,Q_gr,R_gr_td,Qf_gr,Qf_dot_gr,x_ref,th),X0);
toc

z0_star = [0;
          X_grstar(1,1); 
          0; 
          X_grstar(2,1); 
          X_grstar(3,1);
          X_grstar(4,1)];

u_star_td = X_grstar(n_free+1:n_free+length(U0_td),1);  % Optimal Ground input sequence for touchdown
u_star_b = X_grstar(n_free+length(U0_td)+1:end,1);      % Optimal Ground input sequence for brake
U_grstar = [u_star_td; u_star_b];                    % Overall Optimal Ground input sequence

z_ref = z0_star;    % Reference built for the In-Flight optimization

%% Save the result for the new (improved) initial guess
% Xstar = X_grstar;
% save initialguess_gr.mat Xstar

%% Flight Optimization parameters
myoptions   =   myoptimset;
myoptions.ls_beta       = 0.3;        
% myoptions.ls_c          = 0.1;
% myoptions.tolgrad    	=	1e-6;
% myoptions.tolfun    	=	1e-18;
% myoptions.tolX       	=	1e-18;
% normal
myoptions.gradmethod    = 'CD';
myoptions.graddx        = 2^-17;
% user provided
% myoptions.gradmethod    = 'UP';
% myoptions.graddx        = eps^(1/3);
myoptions.nitermax      = 50;
myoptions.Hessmethod    = 'GN';
myoptions.BFGS_gamma  	= 1e-1;
myoptions.GN_sigma      = 1e-3;
myoptions.GN_funF       = @(X_fl)Flight_cost_GN(X_fl,z0_fl,d,Ts,Tend_fl,ds_u,Q_fl,R_fl,Qf_fl,z_ref,th);

%% Self-made initial guess for the Flight optimization
   
U0_fl          = 1*ones(4*N_fl/ds_u,1);
for i = 1:length(U0_fl)
    if mod(i,4) == 1
        U0_fl(i,1) = 0.80;
    elseif mod(i,4) == 2
        U0_fl(i,1) = 0.1;
    elseif mod(i,4) == 3
        U0_fl(i,1) = 0.1;
    elseif mod(i,4) == 0
        U0_fl(i,1) = 0;
    end
end

X0_fl = U0_fl;

%% Improving initial guess for the Flight optimization
load('initialguess_fl.mat');
X0_fl = Xstar;

%% Generating code (mex files) for faster computation (to be ran only once whenever Tend changes)
% codegen Pendulum_cost -args {U,z0,d,Ts,Tend,Q,R,Qf,th}
% codegen Flight_cost_GN -args {U0,z0,Ts,Q,R,Qf,parameters}

%% Running the Flight optimization routine
tic
[X_flstar,fxstar_fl] = myfminunc(@(X_fl)Flight_cost(X_fl,z0_fl,d,Ts,Tend_fl,ds_u,Q_fl,R_fl,Qf_fl,z_ref,th),X0_fl,myoptions);
% [X_flstar,fxstar_fl] = fminunc(@(X_fl)Flight_cost(X_fl,z0_fl,d,Ts,Tend_fl,ds_u,Q_fl,R_fl,Qf_fl,z_ref,th),X0_fl);
toc

U_flstar  = X_flstar;

%% Simulation of Results
Ntot = N_fl+N_gr;
z = zeros(nz,Ntot+1);
zd = zeros(nz,Ntot+1);
z(:,1) = z0_fl;
i = 1;
u_gr = zeros(nu_gr,N_gr);
u_fl = zeros(nu_fl,N_fl);
u = zeros(7,Ntot);


% RK2 simulation
% Flight
u_check = zeros(nu_gr+nu_fl,1);
for ind = 2:N_fl+1
    u_check(nu_gr+1:end,1) = U_flstar(i:i+nu_fl-1,1);
    if mod(ind,ds_u) == 0
        if i < Tend_fl/(Ts*ds_u)*nu_fl-nu_fl
            i = i+nu_fl;
        end
    end
    z_dot = fly2(0, z(:,ind-1), u_check(nu_gr+1:end,1), d, th);
    u(1:nu_fl,ind-1) = u_check(nu_gr+1:end,1);
    zd(:,ind-1) = z_dot;
    z_prime = z(:,ind-1) + Ts/2*z_dot;
    z(:,ind) = z(:,ind-1) + Ts*fly2(0, z_prime, u_check(nu_gr+1:end,1), d, th);
end

% Ground
u_check = zeros(nu_gr+nu_fl,1);

for i = 1:N_td/ds_u_td
    u_check(1:nu_gr,1) = u_star_td((i-1)*nu_gr+1:i*nu_gr,1);
    for ind = N_fl+(i-1)*ds_u_td+2:N_fl+i*ds_u_td+1
        z_dot = ground2(0, z(:,ind-1), u_check(1:nu_gr,1), d, th);
        u(1,ind-1) = u_check(1,1);
        u(2,ind-1) = u_check(2,1);
        u(3,ind-1) = u_check(3,1);
        u(5,ind-1) = u_check(4,1);
        u(6,ind-1) = u_check(5,1);
        u(7,ind-1) = u_check(6,1);
        zd(:,ind-1) = z_dot;
        z_prime = z(:,ind-1) + Ts/2*z_dot;
        z(:,ind) = z(:,ind-1) + Ts*ground2(0, z_prime, u_check(1:nu_gr,1), d, th);
    end
end

for i = 1:N_b/ds_u_b
    u_check(1:nu_gr,1) = u_star_b((i-1)*nu_gr+1:i*nu_gr,1);
    for ind = N_fl+N_td-1+(i-1)*ds_u_b+2:N_fl+N_td+i*ds_u_b+1
        z_dot = ground2(0, z(:,ind-1), u_check(1:nu_gr,1), d, th);
        u(1,ind-1) = u_check(1,1);
        u(2,ind-1) = u_check(2,1);
        u(3,ind-1) = u_check(3,1);
        u(5,ind-1) = u_check(4,1);
        u(6,ind-1) = u_check(5,1);
        u(7,ind-1) = u_check(6,1);
        zd(:,ind) = z_dot;
        z_prime = z(:,ind-1) + Ts/2*z_dot;
        z(:,ind) = z(:,ind-1) + Ts*ground2(0, z_prime, u_check(1:nu_gr,1), d, th);
    end
end

save z 
save zd
save u

%% Plot of the States

time = 0:Ts:Tend;

figure(1)
hold on
plot(time,z(1,:));
plot(time,zd(1,:));
xlabel('Time [s]','Interpreter','latex','FontSize',13);
grid 
legend('$X$','$\dot{X}$','Interpreter','latex','FontSize',13);
title('Longitudinal position and velocity','Interpreter','latex','FontSize',13);
hold off

figure(2)
hold on
plot(time,z(2,:));
plot(time,zd(2,:));
xlabel('Time [s]','Interpreter','latex','FontSize',13);
grid 
legend('$\dot{X}$','$\ddot{X}$','Interpreter','latex','FontSize',13);
title('Longitudinal speed and acceleration','Interpreter','latex','FontSize',13)
hold off

figure(3)
hold on
plot(time,z(3,:));
plot(time,zd(3,:));
xlabel('Time [s]','Interpreter','latex','FontSize',13);
grid 
legend('$Z$','$\dot{Z}$','Interpreter','latex','FontSize',13);
title('Vertical position and velocity','Interpreter','latex','FontSize',13)
hold off

figure(4)
hold on
% plot(time,z(4,:),'DisplayName','speed');
plot(time,zd(4,:),'LineWidth',1.5);
grid 
title('Vertical acceleration','Interpreter','latex','FontSize',13)
xlabel('Time [s]','Interpreter','latex','FontSize',13);
ylabel('$\ddot{Z}$ [$m/s^{2}$]','Interpreter','latex','FontSize',13);
hold off

clear t
size=26;
figure(15)
t=tiledlayout(3,1);
t(1)=nexttile;

plot(t(1),time,180/pi*zd(6,:),"LineWidth",1.5);
grid 
title('Pitch acceleration','Interpreter','latex','FontSize',size)
ylabel('$\ddot{\theta}$ [$deg/s^{2}$]','Interpreter','latex','FontSize',size);
ax = gca; 
ax.XAxis.FontSize = 20; 
ax.YAxis.FontSize = 20;

t(2)=nexttile;
plot(t(2),time,zd(4,:),"LineWidth",1.5);
grid 
title('Vertical acceleration','Interpreter','latex','FontSize',size)
ylabel('$\ddot{Z}$ [$m/s^{2}$]','Interpreter','latex','FontSize',size);
ax = gca; 
ax.XAxis.FontSize = 20; 
ax.YAxis.FontSize = 20;

t(3)=nexttile;
plot(t(3),time,zd(2,:),"LineWidth",1.5);
grid
title('Longitudinal acceleration','Interpreter','latex','FontSize',size)
xlabel('Time [s]','Interpreter','latex','FontSize',size);
ylabel('$\ddot{X}$ [$m/s^{2}$]','Interpreter','latex','FontSize',size);
ax = gca; 
ax.XAxis.FontSize = 20; 
ax.YAxis.FontSize = 20;



figure(5)
hold on
plot(time,180/pi*z(5,:),'LineWidth',1.5);
% plot(time,180/pi*zd(5,:),'r','DisplayName','speed');
grid 
title('Pitch angle','Interpreter','latex','FontSize',13)
xlabel('Time [s]','Interpreter','latex','FontSize',13);
ylabel('$\theta$ [deg]','Interpreter','latex','FontSize',13)
hold off

figure(6)
hold on
% plot(time,180/pi*z(6,:),'DisplayName','speed');
plot(time,180/pi*zd(6,:),'LineWidth',1.2);
grid 
title('Pitch acceleration','Interpreter','latex','FontSize',13)
xlabel('Time [s]','Interpreter','latex','FontSize',13);
ylabel('$\ddot{\theta}$ [$deg/s^{2}$]','Interpreter','latex','FontSize',13);
hold off

figure(7)
plot(z(1,:), z(3,:),'LineWidth',1.5);
grid
title('Trajectory','Interpreter','latex','FontSize',13)
xlabel('X [m]','Interpreter','latex','FontSize',13);
ylabel('Y [m]','Interpreter','latex','FontSize',13);

%% Plot of the Inputs

figure(8)
plot(time(1:end-1),T_max*u(1,:));
grid
title('Thrust','Interpreter','latex','FontSize',13);
xlabel('Time [s]','Interpreter','latex','FontSize',13);
ylabel('$T$ [N]','Interpreter','latex','FontSize',13);

figure(9)
plot(time(1:end-1),u(2,:));
grid
title('Lift Coefficient','Interpreter','latex','FontSize',13);
xlabel('Time [s]','Interpreter','latex','FontSize',13);
ylabel('$u_L$','Interpreter','latex','FontSize',13);

figure(10)
plot(time(1:end-1),u(3,:));
grid
title('Drag Coefficient','Interpreter','latex','FontSize',13);
xlabel('Time [s]','Interpreter','latex','FontSize',13);
ylabel('$u_D$','Interpreter','latex','FontSize',13);

figure(11)
plot(time(1:end-1),theta_max/pi*180*u(4,:));
grid
title('Pitch reference','Interpreter','latex','FontSize',13);
xlabel('Time [s]','Interpreter','latex','FontSize',13);
ylabel('$\theta_in$ [deg]','Interpreter','latex','FontSize',13);

figure(12)
plot(time(1:end-1),Brake_max*u(5,:));
grid
title('Brake Force','Interpreter','latex','FontSize',13);
xlabel('Time [s]','Interpreter','latex','FontSize',13);
ylabel('$F_B$ [N]','Interpreter','latex','FontSize',13);

figure(13)
plot(time(1:end-1),Fa_max*u(6,:));
grid
title('Rear Active suspension Force','Interpreter','latex','FontSize',13);
xlabel('Time [s]','Interpreter','latex','FontSize',13);
ylabel('$F_{Ar}$ [N]','Interpreter','latex','FontSize',13);

figure(14)
plot(time(1:end-1),Fa_max*u(7,:));
grid
title('Front Active suspension Force','Interpreter','latex','FontSize',13);
xlabel('Time [s]','Interpreter','latex','FontSize',13);
ylabel('$F_{Afr}$ [N]','Interpreter','latex','FontSize',13);



[f2_Bellman,zsim2]=Ground_cost(X_grstar,n_free,nu_gr,d,Ts,Tend_td,Tend_b,ds_u_td,ds_u_b,Q_gr,R_gr_td,Qf_gr,Qf_dot_gr,x_ref,th);
f2_Bellman
[f1_Bellman,zsim1]=Flight_cost(X_flstar,z0_fl,d,Ts,Tend_fl,ds_u,Q_fl,R_fl,Qf_fl,z_ref,th);
f1_Bellman
f_Bellman=f1_Bellman+f2_Bellman


%% Animation
% figure(12)
% for i=1/Ts:1:Tend/Ts
%     grid
% %     subplot 311
% %     hold on
% %     plot(z(1,i),z(2,i),'k*', 'linewidth',0.1);
% %     axis([8450 9150 60 85]);
% %     hold off
% %     
% %     subplot 312
% %     hold on
% %     plot(z(1,i),z(4,i),'b*', 'linewidth',0.1);
% %     axis([8450 9150 -10 5]);
% %     hold off
% %     
% %     subplot 313
% %     plot([z(1,i)-Lr*cos(z(5,i)) z(1,i)+Lf*cos(z(5,i))],...
% %         [z(3,i)-Lr*sin(z(5,i)) z(3,i)+Lf*sin(z(5,i))],'-k*', 'linewidth',1);
% %     title(['Time: ' num2str((i-1)*Ts) ' s'])
% %     axis([z(1,i)-20 z(1,i)+25 z(3,i)-10 z(3,i)+10]);
% %     pause(1e-3);
%       plot([z(1,i)-Lr*cos(z(5,i)) z(1,i)+Lf*cos(z(5,i))],...
%         [z(3,i)-Lr*sin(z(5,i)) z(3,i)+Lf*sin(z(5,i))],'-k*', 'linewidth',1);
%       title(['Time: ' num2str((i-1)*Ts) ' s'])
%       axis([z(1,i)-20 z(1,i)+25 z(3,i)-2 z(3,i)+2]);
%       pause(1e-3);
%       
% end

%% barrier tuning for the lift drag input
% ul_max = 1;
% ul_min = 0;
% alpha = 0.01;
% beta = 7e2;
% figure
% i = -0.01:0.0001:1.01;
% plot(i,alpha^2*exp(-2*beta*(ul_max - i)) + alpha^2*exp(-2*beta*(i - ul_min)));

% %% barrier tuning for the thrust
% T_max = 60000;
% T_min = 0;
% 
% alpha = 1;
% beta = 0.025;
% figure
% i = -100:1:60100;
% plot(i,alpha^2*exp(-2*beta*(T_max - i)) + alpha^2*exp(-2*beta*(i - T_min)));
