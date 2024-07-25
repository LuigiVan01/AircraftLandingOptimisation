function [zdot, data_ground] = ground2(t, z, u, d, th)

%% definisci costanti th=[]; kground? where we can put it?
M       = th(1,1);      % mass aircraft
J       = th(2,1);      % inertia
m       = th(3,1);      % wheel mass
k_r     = th(4,1);      % suspension stiffness
k_f     = th(5,1);      % wheel stiffness
c       = th(6,1);      % suspension damping
c_w     = th(7,1);      % wheel damping
S       = th(8,1);      % aircraft surface
L       = th(9,1);      % aircraft length
Lr      = th(10,1);     % rear length
Lf      = th(11,1);     % front length
T_max   = th(12,1);
theta_max = th(13,1);
Brake_max = th(14,1);
Fa_max  = th(15,1);

rho     = 1.225;        % air density [kg/m^3]
g       = 9.81;         % gravitational acceleration [m/s^2]

z_bar = 0;
theta_bar = 0;
M_ver = M+3*m;

A_eq = [  2*k_r    k_f;
        -2*Lr*k_r Lf*k_f];
b_eq = [z_bar*(2*k_r + k_f) + sin(theta_bar)*(Lf*k_f - 2*Lr*k_r) + M_ver*g;
        z_bar*(Lf*k_f - 2*Lr*k_r) - sin(theta_bar)*(Lf^2*k_f + 2*Lr^2*k_r)];
        
Delta = A_eq\b_eq;
Dr = Delta(1,1);
Df = Delta(2,1);

cl0     = 0.63;         % lift coeff. when alpha=0
cl_alpha_max = 3.15;
cl_alpha_var = 0.2;     % changed, the value we set before was 0.8
cd0      = 0.12;         % lift coeff. when alpha=0
gamma    = 0.0435;



%% definisci stati
X          = z(1,1);    % inertial X position (m)
X_dot      = z(2,1);    % body x velocity (m/s)
Z          = z(3,1);    % aircraft height (m)
Z_dot      = z(4,1);    % body z velocity (m/s)
theta      = z(5,1);    % inertial Y position (m)
theta_dot  = z(6,1);    % pitch velocity (rad/s)


%% definisci inputs
uT     =     u(1,1);
uL     =     u(2,1);    % flap opening 
uD     =     u(3,1);    % air-brakes opening
uB     =     u(4,1);    % brake force(N)
uA_r   =     u(5,1);    % rear active suspension force
uA_f   =     u(6,1);    % front active suspension force

T = T_max*uT;
F_brake = Brake_max*uB;
Fa_r = uA_r*Fa_max;
Fa_f = uA_f*Fa_max;

%% Disturbances
v_wind =     d(1,1);    % wind velocity               

%% Equazioni statiche

% without theta_dot effects
v_rel = sqrt((v_wind - X_dot)^2 + (-Z_dot)^2);  
Psi = atan2(-Z_dot, v_wind - X_dot);
Psi_a = atan((-Z_dot)/(v_wind - X_dot));
alpha = theta - Psi_a;

cl = cl0 + (cl_alpha_var + (1-cl_alpha_var)*uL)*cl_alpha_max*alpha;
cd = cd0*(1+uD)+gamma*cl^2;
Fl = 0.5*rho*S*cl*(v_rel)^2;
Fd = 0.5*rho*S*cd*(v_rel)^2;

M_ver = M+3*m;      % Aircraft total mass


Fs_r        = 2*k_r*(Z-Lr*sin(theta)-Dr)+2*Fa_r...
                  +2*c*(Z_dot-Lr*cos(theta)*theta_dot); % rear suspensions total force
Fs_fr       = k_f*(Z+Lf*sin(theta)-Df)+Fa_f...
                  +c*(Z_dot+Lf*cos(theta)*theta_dot);  % front suspension total force
                                                                       
X_ddot      = -(F_brake)/M_ver+Fl/M_ver*cos(Psi-pi/2)...          % longitudinal acceleration
                  +Fd/M_ver*cos(Psi)+T/M_ver*cos(theta);                    
Z_ddot      = -g-Fs_r/M-Fs_fr/M+Fd/M*sin(Psi)+...        % vertical acceleration
                   Fl/M*sin(Psi-pi/2)-T/M*sin(theta);       
theta_ddot  = -Lf/J*Fs_fr+Lr/J*Fs_r;                    % pitch angle acceleration

%% equazioni dinamiche del modello

zdot(1,1)  = X_dot;             % longitudinal velocity
zdot(2,1)  = X_ddot;            % longitudinal acceleration
zdot(3,1)  = Z_dot;             % vertical velocity
zdot(4,1)  = Z_ddot;            % vertical acceleration
zdot(5,1)  = theta_dot;         % pitch angle rate
zdot(6,1)  = theta_ddot;        % pitch angle acceleration

data_ground = [Fl;
               Fd;
               Fs_r;
               Fs_fr];

end
