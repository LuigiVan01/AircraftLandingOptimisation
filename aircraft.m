function [zdot, F] = aircraft(t, z, u, d, th)

%% definisci costanti th=[]; kground? where we can put it?
M   = th(1,1);          % mass aircraft
J   = th(2,1);          % inertia
m   = th(3,1);          % wheel mass
k   = th(4,1);          % suspension stiffness
k_w = th(5,1);          % wheel stiffness
c   = th(6,1);          % suspension damping
c_w = th(7,1);          % wheel damping
S   = th(8,1);          % aircraft surface
L   = th(9,1);          % aircraft length
Lr  = th(10,1);         % rear length
Lf  = th(11,1);         % front length
% Cl0, Cd0, Cm0, deltaCl,...? where to put...
% rho = ;              %air density
g=9.81;

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
z_g_fr      = z(12,1);  % front ground height(m) 


%% definisci input
T      =       u(1,1);   % thrust force(N)
F_r    =       u(2,1);   % brake rear force(N)
F_fr   =       u(3,1);   % brake front force(N)
uL     =       u(4,1);   % flap opening 
uD     =       u(5,1);   % air-brakes opening 
Fa_r   =       u(6,1);   % rear active suspension force
Fa_f   =       u(7,1);   % front active suspension force
wind_v =       d(1,1);   % wind velocity
% come consideriamo il momento?

%% equazioni statiche
%L=
%D=
%M=
Fs_r=2*k*(Z-z_w_r-Lr*sin(theta)+DL0)+2*Fa_r+2*c*(Z_dot-z_w_r_dot-Lr*cos(theta)*theta_dot);   % rear suspensions total force
Fs_fr=k*(Z-z_w_fr+Lf*sin(theta)+DL0)+Fa_f+c*(Z_dot-z_w_f_dot+Lf*cos(theta)*theta_dot);       % front suspension total force
Fw_r=2*k_w*(z_w_r-z_g_r+DeltaS)+2*c_w*(z_w_r_dot-z_g_r_dot);                                 % rear wheel force
Fw_fr=k_w*(z_w_fr-z_g_fr+DeltaS)+c_w*(z_w_fr_dot-z_g_fr_dot);                                % front wheel force
%k_g_r=
%k_g_fr=

%% equazioni dinamiche del modello

zdot(1,1)=X_dot;                                                                    % longitudinal velocity
zdot(2,1)=-(F_fr+2*F_r)/M-L/M*sin(theta)-D/M*cos(theta)+T/M*cos(theta);             % longitudinal acceleration
zdot(3,1)=Z_dot;                                                                    % vertical velocity
zdot(4,1)=-g-Fs_r/M-Fs_fr/M-D/M*sin(theta)+L/M*cos(theta)-T/M*sin(theta);           % vertical acceleration
zdot(5,1)=theta_dot;                                                                % pitch angle rate
zdot(6,1)=Input_Torque/J-Lf/J*Fs_fr+Lr/J*Fs_r;                                      % pitch angle acceleration
zdot(7,1)=z_w_r_dot;                                                                % rear wheel vertical velocity
zdot(8,1)=-g+Fs_r/(2*m)-Fw_r/(2*m);                                                 % rear wheel vertical acceleration
zdot(9,1)=z_w_fr_dot;                                                               % front wheel vertical velocity
zdot(10,1)=-g+Fs_fr/m-Fw_fr/m;                                                      % front wheel vertical acceleration
zdot(11,1)=(k_g_r*(z_g_r+DeltaSg)-k_w*(z_w_r-z_w_r+DeltaS))/c_w+z_w_r_dot;          % rear ground velocity % non ce bisogno del fattore 2
zdot(12,1)=(k_g_fr*(z_g_fr+DeltaSg)-k_w*(z_w_fr-z_w_fr+DeltaS))/c_w+z_w_fr_dot;     % front ground velocity 




zdot = ; % output vector of state
F    = ; % outpuc vector of forces
end
