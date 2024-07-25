clear all
close all

%% Model parameters

M     = 22e3;               % aircraft mass
m     = 130;                % wheel mass
J     = 100e3;              % inertia
k_f   = 6.73e5;             % front stiffness
k_r   = 1.59e4;             % rear stiffness
c     = 4066;               % suspension damping
c_w   = 1.43e5;             % wheel damping
L     = 10;                 % aircraft length
S     = 1/2*10^2;           % aircraft surface
Lf    = 7.76;               % front length
Lr    = 1.94;               % rear length

T_max     = 6e4;            % maximum thrust
theta_max = 10/180*pi;      % maximum pitch
Brake_max = 4e5;            % maximum brake force
Fa_max    = 1e5;            % maximum active suspension force

th = [M J m k_f k_r c c_w S L Lr Lf T_max theta_max Brake_max Fa_max]';
g = 9.81;

d = 0;

%% Control problem parameters
Ts          =   0.01;                        % Sampling time
Tend_fl     =   15;                          % Time horizon
N_fl        =   Tend_fl/Ts;                  % Prediction steps

% Downsampling the inputs
ds_T        = 50;
ds_L        = 50;
ds_D        = 100;
ds_theta    = 25;
ds_u_fl = [ds_T, ds_L, ds_D, ds_theta]';     

nz          =   6;
nu_fl       =   4;

%% Landing strip length

x_ref    =   [ 1000   ;          % lower bound
               1600 ];           % upper bound

%% Initial condition for the Flight operation
hor_speed_in = 100;
height_in = 50;
vert_speed_in = -5;
pitch_in = 0;
pitch_dot_in = 0;

z0_fl = [0; hor_speed_in; height_in; vert_speed_in; pitch_in; pitch_dot_in];

%% Reference for the final landing target
load('ref_for_flight.mat');            % Results of the Ground optimization
z_ref = z0_ref;

%% Optimization parameters 

zdd_w       =   1;                                  % weight on vertical acceleration
thdd_w      =   1*180/pi;                           % weight on pitch acceleration
xdd_w       =   1;                                  % weight on longitudinal acceleration
Q_fl        =   diag([0;xdd_w;0;zdd_w;0;thdd_w]);   % state weighting matrix
R_fl        =   1;                                  % Input weight: Thrust && "flare" penalization

%% Initial Guess

T0 = 0.6*ones(N_fl/ds_u_fl(1,1),1);
L0 = 1*ones(N_fl/ds_u_fl(2,1),1);
D0 = 1*ones(N_fl/ds_u_fl(3,1),1);
th0 = 0*ones(N_fl/ds_u_fl(4,1),1);

U0 = [T0; L0; D0; th0];
X0 = U0;

%% Build your initial guess
% Use the function new_prova_fl and comment all the following "optimization
% parts" of this script and just plot the results. 
% To modify the initial guess try to act on the "if clauses" regarding the inputs into the
% simulation code of the new_prova_fl script


% [z_sim, z_dot, T0, L0, D0, th0] = ...
%     new_prova_fl(X0,z0_fl,nu_fl,nz,d,Ts,Tend_fl,ds_u_fl,Q_fl,R_fl,z_ref,th);
% 
% z = zeros(nz, N_fl+1);
% z_d = zeros(nz, N_fl+1);
% 
% for ind = 1:N_fl+1
%     z(:,ind) = z_sim((ind-1)*nz+1:ind*nz);
%     z_d(:,ind) = z_dot((ind-1)*nz+1:ind*nz);
% end
% 
% U0 = [T0; L0; D0; th0];
% X0 = U0;

%% Main parameters for the Flight constraints

% Rate limiters: [T; L; D; theta_ref]
rate_bound = [0.4, 0.5, 0.5, 10/20]';           % rate bound for each input

% Upper and Lower bounds: [T; L; D; theta_ref]
ub_input = ones(nu_fl,1);                       % upper bounds for the input sequence
lb_input = [0;0;0;-1];                          % lower bounds for the input sequence

%% Linear equality constraint parameters

A  = [];
b  = [];

%% Linear inequality constraint parameters

% Rate Limiter for each input
C_rate_lim = zeros(2*length(U0),length(X0));    % initialization
C_bound = zeros(2*length(U0),length(X0));       % initialization
d_rate_lim = zeros(2*length(U0),1);             % initialization
d_bound = zeros(2*length(U0),1);                % initialization

Nu_fl = zeros(nu_fl,1);                         % number of samples for each input

l = 0;                                          % ending of each input block in the matrix  
fin = 0;                                        % start of each input block in the matrix
for i = 1:nu_fl
    Nu_fl(i,1) = N_fl/ds_u_fl(i,1);                
    l = l + 2*Nu_fl(i,1);
    
    % Rate limiter
    C_rate_lim(fin+1:l,fin/2+1:l/2) = ...
        gen_rate_lim_mat(Nu_fl(i,1));
    d_rate_lim(fin+1:l,1) = rate_bound(i,1)*ones(2*Nu_fl(i,1),1);
    
    % Upper and lower bound for the input
    C_bound(fin+1:l,fin/2+1:l/2) = ...
        [-eye(Nu_fl(i,1)); eye(Nu_fl(i,1))];
    d_bound(fin+1:l,1) = [-ub_input(i,1)*ones(Nu_fl(i,1),1); lb_input(i,1)*ones(Nu_fl(i,1),1)];
    
    fin = fin + 2*Nu_fl(i,1);
end

% Final matrix with rate limiters
C = [C_rate_lim; C_bound];
d = [-d_rate_lim; d_bound];


%% Nonlinear constraints specifications

% Equality constraints: to enforce the macthing between the final state of
% the flight and the optimal initail condition of the ground.

% Inequality constraints: to forbid the aircraft to have heigth smaller
% than zero. To speed-up the optimization routines, just the final moments
% of the flight are taken into account.


p = nz-1;         % number of NL equality constraints;
q = N_fl/10;      % number of NL inequality contraints;

%% Gauss-Newton
% Initialize solver options
myoptions               =   myoptimset;
myoptions.Hessmethod  	=	'GN';
myoptions.GN_funF       =	@(X_fl)new_Flight_cost_GN(X_fl,z0_fl,nu_fl,nz,d,Ts,Tend_fl,ds_u_fl,Q_fl,R_fl,z_ref,x_ref, th);
myoptions.gradmethod  	=	'CD';
myoptions.graddx        =	2^-17;
myoptions.tolgrad    	=	1e-6;
myoptions.ls_nitermax   =	20;
myoptions.nitermax      =	50;
myoptions.xsequence     =	'on';


tic
[xstar,fxstar,niter,exitflag,xsequence] = ...
    myfmincon(@(X_fl)new_Flight_cost(X_fl,z0_fl,nu_fl,nz,d,Ts,Tend_fl,ds_u_fl,Q_fl,R_fl,z_ref,th),...
    X0,A,b,C,d,p,q,myoptions);
toc

%% Results

% Simulation with the optimal input sequence
[z_sim, z_dot, Tstar, Lstar, Dstar, thstar, data_fly] = ...
    new_sim_fl(xstar,z0_fl,nu_fl,nz,d,Ts,Tend_fl,ds_u_fl,Q_fl,R_fl,z_ref,th);

z = zeros(nz, N_fl+1);
z_d = zeros(nz, N_fl+1);

for ind = 1:N_fl+1
    z(:,ind) = z_sim((ind-1)*nz+1:ind*nz);
    z_d(:,ind) = z_dot((ind-1)*nz+1:ind*nz);
end

% Save the results
flight.sim = [z; z_d];
flight.T = Tstar;
flight.L = Lstar;
flight.D = Dstar;
flight.Theta_ref = thstar;
flight.Fl = data_fly(1,:);
flight.Fd = data_fly(2,:);
flight.alpha = data_fly(3,:);
flight.v_rel = data_fly(4,:);

save Result_Flight.mat flight

%% Plot of the States
time = 0:Ts:Tend_fl;

figure(1)
hold on
plot(time,z(1,:),'b','DisplayName','position');
plot(time,z_d(1,:),'r','DisplayName','speed');
grid
legend
title('Horizontal position overview',"Interpreter","Latex")
xlabel('Time',"Interpreter","Latex");
hold off

figure(2)
hold on
plot(time,z(2,:),'b','DisplayName','speed');
plot(time,z_d(2,:),'r','DisplayName','acceleration');
grid 
legend
title('Horizontal speed overview',"Interpreter","Latex")
xlabel('Time',"Interpreter","Latex");
hold off

figure(3)
hold on
plot(time,z(3,:),'b','DisplayName','position');
plot(time,z_d(3,:),'r','DisplayName','speed');
grid 
legend
title('Vertical position overview',"Interpreter","Latex")
xlabel('Time',"Interpreter","Latex");
hold off

figure(4)
hold on
plot(time,z(4,:),'b','DisplayName','speed');
plot(time,z_d(4,:),'r','DisplayName','acceleration');
grid
legend
title('Vertical acceleration overview',"Interpreter","Latex")
xlabel('Time',"Interpreter","Latex");
hold off

figure(5)
hold on
plot(time,180/pi*z(5,:),'b','DisplayName','position');
plot(time,180/pi*z_d(5,:),'r','DisplayName','speed');
grid 
legend
title('Pitch overview',"Interpreter","Latex")
xlabel('Time',"Interpreter","Latex");
hold off

figure(6)
hold on
plot(time,180/pi*z(6,:),'b','DisplayName','speed');
plot(time,180/pi*z_d(6,:),'r','DisplayName','acceleration');
grid 
legend
title('Pitch acceleration overview',"Interpreter","Latex")
xlabel('Time',"Interpreter","Latex");
hold off

figure(7)
plot(z(1,:), z(3,:),'k');
grid
title('Trajectory',"Interpreter","Latex")

%% Plots of the inputs
fin = 0;
for i = 1:length(Tstar)
    T_real(fin+1:fin+ds_u_fl(1,1),1) = repmat(Tstar(i,1),ds_u_fl(1,1),1);
    fin = fin + ds_u_fl(1,1);
end

fin = 0;
for i = 1:length(Lstar)
    L_real(fin+1:fin+ds_u_fl(2,1),1) = repmat(Lstar(i,1),ds_u_fl(2,1),1);
    fin = fin + ds_u_fl(2,1);
end

fin = 0;
for i = 1:length(Dstar)
    D_real(fin+1:fin+ds_u_fl(3,1),1) = repmat(Dstar(i,1),ds_u_fl(3,1),1);
    fin = fin + ds_u_fl(3,1);
end

fin = 0;
for i = 1:length(thstar)
    th_real(fin+1:fin+ds_u_fl(4,1),1) = repmat(thstar(i,1),ds_u_fl(4,1),1);
    fin = fin + ds_u_fl(4,1);
end

figure(8)
plot(T_max*T_real);
title('Thrust',"Interpreter","Latex");
xlabel('Time',"Interpreter","Latex");

figure(9)
plot(L_real);
title('Lift Coefficient',"Interpreter","Latex");
xlabel('Time',"Interpreter","Latex");

figure(10)
plot(D_real);
title('Drag Coefficient',"Interpreter","Latex");
xlabel('Time',"Interpreter","Latex");

figure(11)
plot(theta_max/pi*180*th_real);
title('Pitch reference',"Interpreter","Latex");
xlabel('Time',"Interpreter","Latex");



     



