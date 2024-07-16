clear all
close all
clc

%% Model parameters

M = 22e3;
m = 130;
J = 100e3;
k_f = 6.73e5;
k_r = 1.59e4;            % modified: different (smaller) value for the front stiffness
% k_r = 6.73e5;
c = 4066;
c_w = 1.43e5;
L = 10;
S = 1/2*10^2;
Lf = 7.76;
Lr = 1.94;
T_max = 6e4;
theta_max = 10/180*pi;   % modified
Brake_max = 4e5;         
Fa_max = 1e5;

th = [M J m k_f k_r c c_w S L Lr Lf T_max theta_max Brake_max Fa_max]';
g = 9.81;

d = 0;

%% Control problem parameters
Ts          =   0.01;                        % Sampling time
Tend_gr     =   15;                          % Time horizon
N_gr        =   Tend_gr/Ts;                  % Prediction steps

ds_T        = 500;
ds_L        = 150;
ds_D        = 150;
ds_B        = 75;
ds_Far      = 15;
ds_Faf      = 15;
ds_u_gr = [ds_T, ds_L, ds_D, ds_B, ds_Far, ds_Faf]';     % Downsampling the inputs

nz          =   6;
nu_gr       =   6;

%% Landing strip length
x_ref       =   [ 1700  ;       % upper bound
                    0  ];       % lower bound

%% Optimization parameters 
zdd_w       =   1;                                  % weight on vertical acceleration
thdd_w      =   1*180/pi;                           % weight on pitch acceleration
Q_dot_gr    =   diag([0;1;0;zdd_w;0;thdd_w]);       % acceleration weighting matrix
Q_gr        =   diag([0;0;sqrt(1000);0;0;0]);       % state weighting matrix
R_gr        =   0;                                  % Active suspension weight
n_free      =   4;                                  % Degree of freedom in the initial condition optimization

%% Initial Guess

z0_free = [95; -0.5; 3*pi/180; -3*pi/180];          % [hor_speed; vert_speed; pitch; pitch_speed]
T0 = 0*ones(N_gr/ds_u_gr(1,1),1);
L0 = 0*ones(N_gr/ds_u_gr(2,1),1);
D0 = 0*ones(N_gr/ds_u_gr(3,1),1);
B0 = 0*ones(N_gr/ds_u_gr(4,1),1);
Far0 = zeros(N_gr/ds_u_gr(5,1),1);
Faf0 = zeros(N_gr/ds_u_gr(6,1),1);

U0 = [T0; L0; D0; B0; Far0; Faf0];
X0 = [z0_free; U0];

%% Main parameters for Ground constraints

% Rate limiter: [T; L; D; B; Far; Faf]
rate_bound = [0.1, 0.5, 0.5, 0.3, 0.3, 0.3]';   % rate bound for each input
rate_start = [0.1, 0.5, 0.5, 0.3, 1, 1]';       % rate bound for each input at start

% Upper and Lower bounds for the inputs : [T; L; D; B; Far; Faf]
ub_input = [1;1;1;1;1;1];                       % upper bounds for the input sequence
lb_input = [0;0;0;0;-1;-1];                     % lower bounds for the input sequence

% Upper and lower bound for the initial condition: [hor_speed; vert_speed; pitch; pitch_speed]
ub = [95; -0.3; 10*pi/180; 5*pi/180];
lb = [90; -5; 2*pi/180; -5*pi/180];

%% Linear equality constraint definition
A = [];
b = [];

%% Linear inequality constraint definition
% C = zeros(2*length(U0)+2*length(U0)+4,length(X0));
% d = zeros(2*length(U0)+2*length(U0)+2*length(z0),1);

% Rate Limiter for each input
C_rate_lim  = zeros(2*length(U0),length(X0));       % initialization
C_bound     = zeros(2*length(U0),length(X0));       % initialization
C_start     = zeros(2*nu_gr,length(X0));            % initialization
d_rate_lim  = zeros(2*length(U0),1);                % initialization
d_bound     = zeros(2*length(U0),1);                % initialization
d_start     = zeros(2*nu_gr,1);                     % initialization

Nu_gr = zeros(nu_gr,1);                         % Vector containing the number of samples for each input

l = 0;                                          % ending of each input group in the matrix  
fin = 0;                                        % start of each input group in the matrix

% One iteration for each input
for i = 1:nu_gr
    Nu_gr(i,1) = N_gr/ds_u_gr(i,1);             % Compute the number of samples for the input 
    l = l + 2*Nu_gr(i,1);                       % Ending index
    
    %% Rate limiter

    % Compute the matrix representing the rate limiter for the specific input group
    C_rate_lim(fin+1:l,length(z0_free)+fin/2+1:length(z0_free)+l/2) = ...
        gen_rate_lim_mat(Nu_gr(i,1)); 

    % Compute the matrix representing the known term of the rate limiter 
    d_rate_lim(fin+1:l,1) = rate_bound(i,1)*ones(2*Nu_gr(i,1),1); 

    %% Lower and upper bounds on the inputs

    % Upper and lower bound for the input
    C_bound(fin+1:l,length(z0_free)+fin/2+1:length(z0_free)+l/2) = ...
        [-eye(Nu_gr(i,1)); eye(Nu_gr(i,1))];
    d_bound(fin+1:l,1) = [-ub_input(i,1)*ones(Nu_gr(i,1),1); lb_input(i,1)*ones(Nu_gr(i,1),1)];
    
    %% Rate limiter on the first 

    % starter rate limiter bound
    C_start(i,length(z0_free)+fin/2+1) = 1;
    C_start(i+nu_gr,length(z0_free)+fin/2+1) = -1;
    d_start(i,1) = rate_start(i,1);
    d_start(i+nu_gr,1) = rate_start(i,1);
    
    fin = fin + 2*Nu_gr(i,1);
end

% Matrix representing the upper and lower bound for the initial conditions

C_in_cond = [-eye(n_free), zeros(n_free,length(U0));
              eye(n_free), zeros(n_free,length(U0))]; 
         
% Final matrix with rate limiter

C = [C_rate_lim; 
     C_start; 
     C_bound; 
     C_in_cond];

d = [-d_rate_lim; 
     -d_start; 
      d_bound; 
      -ub; 
       lb;];

% Final matrix without rate limiter
% C = [C_bound; C_in_cond];
% d = [d_bound; -ub; lb;];

%% Nonlinear constraints specifications
p = nz+2;      % number of NL equality constraints;
q = 2;         % number of NL inequality contraints;

%% Gauss-Newton
% Initialize solver options
myoptions               =   myoptimset;
myoptions.Hessmethod  	=	'GN';
myoptions.GN_funF       =	@(X_gr)new_Ground_cost_GN(X_gr,n_free,nu_gr,nz,d,Ts,Tend_gr,ds_u_gr,Q_dot_gr,Q_gr,R_gr,x_ref,th);
myoptions.gradmethod  	=	'CD';
myoptions.graddx        =	2^-17;
myoptions.tolgrad    	=	1e-6;
myoptions.tolfun    	=	1e-12;
myoptions.ls_nitermax   =	20;
myoptions.nitermax      =	50;
myoptions.xsequence     =	'on';



%myoptions.QPoptions     = optimoptions('quadprog','Algorithm','active-set');

% BFGS options
% myoptions.ls_tkmax      =	1;          
% myoptions.ls_beta       =	0.5;
% myoptions.ls_c          =	0.1;

tic
[xstar,fxstar,niter,exitflag,xsequence] = ...
    myfmincon(@(X_gr)new_Ground_cost(X_gr,n_free,nu_gr,nz,d,Ts,Tend_gr,ds_u_gr,Q_dot_gr,Q_gr,R_gr,x_ref,th),...
    X0,A,b,C,d,p,q,myoptions);
toc

% Build the target for the Flight optimization 
z0_star = xstar(1:n_free,1);
z0_ref = [0; z0_star(1,1); 0; z0_star(2,1); z0_star(3,1); z0_star(4,1)];
save ref_for_flight.mat z0_ref

U_star = xstar(n_free+1:end,1);

%% Debug lin constr (ignore this section)
% Code to check the feasibility of the linear constrained problem

% A_f = C_rate_lim;
% b_f = d_rate_lim;
% Aeq = A;
% beq = b;
% 
% up_b = [ub; zeros(length(U0),1)];
% low_b = [lb; zeros(length(U0),1)];
% 
% s = n_free;
% for i=1:nu_gr
%     up_b(s+1:s+Nu_gr(i,1)) = ub_input(i,1)*ones(Nu_gr(i,1),1);
%     low_b(s+1:s+Nu_gr(i,1)) = lb_input(i,1)*ones(Nu_gr(i,1),1);
%     s = s + Nu_gr(i,1);
% end
% 
% f = zeros(size(X0)); % assumes x0 is the initial point
% Xnew = linprog(f,A_f,b_f,Aeq,beq,low_b,up_b);

%% Debug nonlincon (ignore this section)
% code to find a feasible solution of the full constrained problem

% X0_n = (up_b + low_b)/3;
% options = optimoptions("fmincon","Algorithm","interior-point","Display" ...
%     ,"iter","HessianApproximation","bfgs","ScaleProblem",true,"MAxFunctionEvaluations",5e3);
% 
% nonlcon = @(X_gr)non_linconstr(X_gr,n_free,nu_gr,nz,d,Ts,Tend_gr,ds_u_gr,Q_dot_gr,Q_gr,R_gr,x_ref,th);
% [xstar,fxstar,niter,exitflag,xsequence] = ...
%     fmincon(@(X_gr)0,...
%     X0_n,A_f,b_f,Aeq,beq,low_b,up_b,nonlcon,options);
% X0 = xstar;
%% fmincon solution

% up_b = [ub; zeros(length(U0),1)];
% low_b = [lb; zeros(length(U0),1)];
% 
% s = n_free;
% for i=1:nu_gr
%     up_b(s+1:s+Nu_gr(i,1)) = ub_input(i,1)*ones(Nu_gr(i,1),1);
%     low_b(s+1:s+Nu_gr(i,1)) = lb_input(i,1)*ones(Nu_gr(i,1),1);
%     s = s + Nu_gr(i,1);
% end
% 
% A_f = C_rate_lim;
% b_f = d_rate_lim;
% Aeq = A;
% beq = b;
% 
% nonlcon = @(X_gr)non_linconstr(X_gr,n_free,nu_gr,nz,d,Ts,Tend_gr,ds_u_gr,Q_dot_gr,Q_gr,R_gr,x_ref,th);
% 
% options = optimoptions("fmincon","Algorithm","sqp","Display" ...
%     ,"iter","HessianApproximation","bfgs","ScaleProblem",true,"MAxFunctionEvaluations",5e3);
% 
% [xstar,fxstar,niter,exitflag,xsequence] = ...
%     fmincon(@(X_gr)new_Ground_cost_fmincon(X_gr,n_free,nu_gr,nz,d,Ts,Tend_gr,ds_u_gr,Q_dot_gr,Q_gr,R_gr,x_ref,th),...
%     X0,A_f,b_f,Aeq,beq,low_b,up_b,nonlcon,options);
% 
% z0_star = xstar(1:n_free,1);
% U_star = xstar(n_free+1:end,1);

%% Results

% Simulation of the optimal input sequence
[z_sim, z_dot, Tstar, Lstar, Dstar, Bstar, Farstar, Fafstar, data_ground] = ...
    new_sim(xstar,n_free,nu_gr,nz,d,Ts,Tend_gr,ds_u_gr,Q_dot_gr,R_gr,th);

z = zeros(nz, N_gr+1);
z_d = zeros(nz, N_gr+1);

for ind = 1:N_gr+1
    z(:,ind) = z_sim((ind-1)*nz+1:ind*nz);
    z_d(:,ind) = z_dot((ind-1)*nz+1:ind*nz);
end

% Save the results
ground.sim = [z; z_d];
ground.T = Tstar;
ground.L = Lstar;
ground.D = Dstar;
ground.B = Bstar;
ground.Far = Farstar;
ground.Faf = Fafstar;
ground.Fl = data_ground(1,:);
ground.Fd = data_ground(2,:);
ground.Fs_r = data_ground(3,:);
ground.Fs_fr = data_ground(4,:);

save Result_Ground.mat ground

%% Plot of the States
time = 0:Ts:Tend_gr;

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
    T_real(fin+1:fin+ds_u_gr(1,1),1) = repmat(Tstar(i,1),ds_u_gr(1,1),1);
    fin = fin + ds_u_gr(1,1);
end

fin = 0;
for i = 1:length(Lstar)
    L_real(fin+1:fin+ds_u_gr(2,1),1) = repmat(Lstar(i,1),ds_u_gr(2,1),1);
    fin = fin + ds_u_gr(2,1);
end

fin = 0;
for i = 1:length(Dstar)
    D_real(fin+1:fin+ds_u_gr(3,1),1) = repmat(Dstar(i,1),ds_u_gr(3,1),1);
    fin = fin + ds_u_gr(3,1);
end

fin = 0;
for i = 1:length(Bstar)
    B_real(fin+1:fin+ds_u_gr(4,1),1) = repmat(Bstar(i,1),ds_u_gr(4,1),1);
    fin = fin + ds_u_gr(4,1);
end

fin = 0;
for i = 1:length(Farstar) 
    Far_real(fin+1:fin+ds_u_gr(5,1),1) = repmat(Farstar(i,1),ds_u_gr(5,1),1);
    fin = fin + ds_u_gr(5,1);
end

fin = 0;
for i = 1:length(Fafstar) 
    Faf_real(fin+1:fin+ds_u_gr(6,1),1) = repmat(Fafstar(i,1),ds_u_gr(6,1),1);
    fin = fin + ds_u_gr(6,1);
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

figure(12)
plot(Brake_max*B_real);
title('Brake Force',"Interpreter","Latex");
xlabel('Time',"Interpreter","Latex");

figure(13)
plot(Fa_max*Faf_real);
title('Rear Active suspension Force',"Interpreter","Latex");
xlabel('Time',"Interpreter","Latex");

figure(14)
plot(Fa_max*Faf_real);
title('Front Active suspension Force',"Interpreter","Latex");
xlabel('Time',"Interpreter","Latex");
