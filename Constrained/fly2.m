function [zdot, data] = fly2(t, z, u, d, th)

%% Parameters
M       = th(1,1);      % mass aircraft
J       = th(2,1);      % inertia
m       = th(3,1);      % wheel mass
k       = th(4,1);      % suspension stiffness
k_w     = th(5,1);      % wheel stiffness
c       = th(6,1);      % suspension damping
c_w     = th(7,1);      % wheel damping
S       = th(8,1);      % aircraft surface
L       = th(9,1);      % aircraft length
Lr      = th(10,1);     % rear length
Lf      = th(11,1);     % front length 
T_max   = th(12,1);
theta_max = th(13,1);

rho     = 1.225;        % air density [kg/m^3]
g       = 9.81;         % gravitational acceleration [m/s^2]


cl0     = 0.63;         % lift coeff. when alpha=0
cl_alpha_max = 3.15;
cl_alpha_var = 0.2;     % changed, the value we set before was 0.8
cd0      = 0.12;         % lift coeff. when alpha=0
gamma    = 0.0435;

M_ver=M+3*m;

%% definisci stati
X          = z(1,1);    % inertial X position (m)
X_dot      = z(2,1);    % body x velocity (m/s)
Z          = z(3,1);    % aircraft height (m)
Z_dot      = z(4,1);    % body z velocity (m/s)
theta      = z(5,1);    % inertial Y position (m)
theta_dot  = z(6,1);    % pitch velocity (rad/s)


%% definisci inputs
uT          =     u(1,1);     % thrust force(N)
uL          =     u(2,1);      % flap opening 
uD          =     u(3,1);      % air-brakes opening 
u_theta_in =   u(4,1);
T = uT*T_max;
theta_in = u_theta_in*theta_max;

%% Disturbances
v_wind =     d(1,1);    % wind velocity               

%% Equazioni statiche

% code for attack angle

Psi = atan2(-Z_dot, v_wind - X_dot);
Psi_a = atan((-Z_dot)/(v_wind - X_dot));
alpha = theta - Psi_a;
v_rel = sqrt((v_wind - X_dot)^2 + ...
    (-Z_dot)^2);


cl = cl0 + (cl_alpha_var + (1-cl_alpha_var)*uL)*cl_alpha_max*alpha;
cd = cd0*(1+uD)+gamma*cl^2;
Fl = 0.5*rho*S*cl*(v_rel)^2;
Fd = 0.5*rho*S*cd*(v_rel)^2;

X_ddot = 1/M_ver*(T*cos(theta) + Fl*cos(Psi - pi/2) + Fd*cos(Psi));
Z_ddot = -g + 1/M_ver*(T*sin(theta) + Fl*sin(Psi - pi/2) + Fd*sin(Psi));
xi = 0.9;
wn = 1.5;
theta_ddot  = -2*wn*xi*theta_dot-theta*wn^2 + wn^2*theta_in;  


%% Equazioni dinamiche del modello

zdot(1,1)  = X_dot;             % longitudinal velocity
zdot(2,1)  = X_ddot;            % longitudinal acceleration
zdot(3,1)  = Z_dot;             % vertical velocity
zdot(4,1)  = Z_ddot;            % vertical acceleration
zdot(5,1)  = theta_dot;         % pitch angle rate
zdot(6,1)  = theta_ddot;        % pitch angle acceleration

data = [Fl;
        Fd;
        alpha;
        v_rel];

end
