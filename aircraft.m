function [zdot, F] = aircraft(t, z, u, d, th)

%% definisci costanti th=[]; kground? where we can put it?
M   = th(1,1);          % mass aircraft
J   = th(1,2);          % inertia
m   = th(1,3);          % wheel mass
k   = th(1,4);          % suspension stiffness
k_w = th(1,5);          % wheel stiffness
c   = th(1,6);          % suspension damping
c_w = th(1,7);          % wheel damping
S   = th(1,8);          % aircraft surface
L   = th(1,9);          % aircraft length
Lr  = th(1,10);         % rear length
Lf  = th(1,11);         % front length
% Cl0, Cd0, Cm0, deltaCl,...? where to put...
% rho = ;   %air density

%% definisci stati
X          = z(1,1);    % inertial X position (m)
X_dot      = z(2,1);    % body x velocity (m/s)
Z          = z(3,1);    % aircraft height (m)
Z_dot      = z(4,1);    % body z velocity (m/s)
theta      = z(5,1);    % inertial Y position (m)
theta_dot  = z(6,1);    % pitch velocity (rad/s)
z_w_r      = z(7,1);    % rear wheel height (m)
z_w_r_dot  = z(8,1);    % rear wheel velocity (m)
z_w_fr     = z(9,1);    % front wheel height (m)
z_w_fr_dot = z(10,1);   % front wheel velocity(m/s)
z_g_r      = z(11,1);   % rear ground height(m)
z_g_r_dot  = z(12,1);   % rear ground velocity(m/s)
z_g_f      = z(13,1);   % front ground height(m)
z_g_f_dot  = z(14,1);   % front ground velocity(m/s)


%% definisci input
T      =       u(1,1);     % thrust force(N)
F_r    =       u(2,1);      % brake rear force(N)
F_fr   =       u(3,1);     % brake front force(N)
uL     =       u(4,1);      % flap opening 
uD     =       u(5,1);      % air-brakes opening 
wind_v =       d(1,1);     % wind velocity

%% equazioni statiche

%% equazioni dinamiche del modello
zdot = ; % output vector of state
F    = ; % outpuc vector of fores
end
