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
Brake_max = 1e4;
Fa_max = 1e5;
% aero_c = 7.76; 

th = [M J m k k_w c c_w S L Lr Lf T_max theta_max Brake_max Fa_max]';
g = 9.81;

d = 0;

%% Control problem parameters
Ts          =   0.01;                            % Sampling time
Tend        =   12;
ds_u        =   80;
N           =   Tend/Ts;                         % Prediction steps
h_start     =   50;
z0          =   [0; 95; h_start; -5; 0; 0];      % Initial state
zdd_w  =   1;
thdd_w =   1*180/pi;
Q           =   diag([0;0;0;zdd_w;0;thdd_w]);    % Tracking error weight
Qf          =   diag([0;1e3;1e3;1e3;1e3;1e3]);                  % Terminal weight
R           =   zeros(4*N/ds_u);                  % Input weight
z_ref  =   [0; 100; 0; -0.3; 0.02; -0.01]; % NON UTILIZZATO
for i = 1:length(R)
    if mod(i,4) == 1
        R(i,i) = 0;
    elseif mod(i,4) == 2
        R(i,i) = 0;
    elseif mod(i,4) == 3
        R(i,i) = 0;
    elseif mod(i,4) == 0
        R(i,i) = 0;
    end
end

%% Optimization parameters
myoptions   =   myoptimset;
myoptions.ls_beta       = 0.3;        
% myoptions.ls_c          = 0.1;
% myoptions.tolgrad    	=	1e-6;
% myoptions.tolfun    	=	1e-9;
% myoptions.tolX       	=	1e-12;
% normal
myoptions.gradmethod    = 'CD';
myoptions.graddx        = 2^-17;
% user provided
% myoptions.gradmethod    = 'UP';
% myoptions.graddx        = eps^(1/3);
myoptions.nitermax      = 25;
myoptions.Hessmethod    = 'GN';
myoptions.BFGS_gamma  	= 1e-1;
myoptions.GN_sigma      = 1e-3;
myoptions.GN_funF       = @(X)Flight_cost_GN(X,z0,d,Ts,Tend,ds_u,Q,R,Qf,z_ref,th);

%% Initial guess

U0          = 1*ones(4*N/ds_u,1);
for i = 1:length(U0)
    if mod(i,4) == 1
        U0(i,1) = 0.80;
    elseif mod(i,4) == 2
        U0(i,1) = 0.9;
    elseif mod(i,4) == 3
        U0(i,1) = 0.1;
    elseif mod(i,4) == 0
        U0(i,1) = 7/180*pi;
    end
end

load('initialguess_fl.mat');
X0 = Xstar;

%% Generating code (mex files) for faster computation (to be ran only once whenever Tend changes)
% codegen Pendulum_cost -args {U,z0,d,Ts,Tend,Q,R,Qf,th}
% codegen Flight_cost_GN -args {U0,z0,Ts,Q,R,Qf,parameters}

%% Running the optimization routine
tic
[Xstar,fxstar] = myfminunc(@(X)Flight_cost(X,z0,d,Ts,Tend,ds_u,Q,R,Qf,z_ref,th),X0,myoptions);
% [Xstar,fxstar] = fminunc(@(X)Flight_cost(X,z0,d,Ts,Tend,ds_u,Q,R,Qf,z_ref,th),X0);
toc
% [Ustar,fxstar,k,exitflag,xsequence] = myfminunc(@(U)Pendulum_cost_mex(U,z0,Ts,Q,R,Qf,parameters),U0,myoptions); % if using compiled version

%% Save the result for the new (improved) initial guess
% save initialguess_fl.mat Xstar

%% Results
close all
nz = 6;
z = zeros(6,N+1);
zd = zeros(6,N+1);
z(:,1) = z0;
i = 1;
time = 1:N+1;
u = zeros(4,N);

F = Flight_cost_GN(Xstar,z0,d,Ts,Tend,ds_u,Q,R,Qf,z_ref,th);
F_state = F(1:nz*N,1);
F_final_state = F(nz*N+1:nz*(N+1),1);
F_barrier = F(nz*(N+2)+1:end,1);

[~,z_sim] = Flight_cost_GN(Xstar,z0,d,Ts,Tend,ds_u,Q,R,Qf,z_ref,th);
zsim = zeros(nz,N+1);
for ind = 1:N+1
    zsim(:,ind) = z_sim((ind-1)*nz+1:ind*nz,1);
end

% RK2 simulation
for ind = 2:N+1
    if mod(ind,ds_u) == 0
        if i < Tend/(Ts*ds_u)*4-4
            i = i+4;
        end
    end
    z_dot = fly2(0, z(:,ind-1), Xstar(i:i+3,1), d, th);
    u(:,ind-1) = Xstar(i:i+3,1);
    zd(:,ind) = z_dot;
    z_prime = z(:,ind-1) + Ts/2*z_dot;
    z(:,ind) = z(:,ind-1) + Ts*fly2(0, z_prime, Xstar(i:i+3,1), d, th);
end

for j=1:6
    figure(j)
    hold on
    if j == 5 || j == 6
        plot(time,180/pi*z(j,:));
        plot(time,180/pi*zsim(j,:), 'b*');
        plot(time,180/pi*zd(j,:));
    else
        plot(time,z(j,:));
        plot(time,zsim(j,:), 'b*');
        plot(time,zd(j,:));
    end
    hold off
end
for j=1:4
    figure(6+j)
    if j == 4
        plot(20*u(j,:));
    else
        plot(u(j,:));
    end
end
figure(11)
hold on
plot(z(1,:), z(3,:),'g');
plot(zsim(1,:), zsim(3,:), 'b*');
hold off

[f1,z_sim1]=Flight_cost(Xstar,z0,d,Ts,Tend,ds_u,Q,R,Qf,z_ref,th);
f1

%% barrier tuning

alpha = 0.00001;
beta = 3e2;
figure
i = -0.05:0.001:0.05;
plot(i,alpha^2*exp(-2*beta*i) + alpha^2*exp(2*beta*i));

% %% Animation
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
%       axis([z(1,i)-20 z(1,i)+25 z(3,i)-10 z(3,i)+10]);
%       pause(1e-3);
%       
% end