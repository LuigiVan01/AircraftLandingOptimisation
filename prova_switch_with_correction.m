clear all                                                                      
close all
clc

%% Parameters
M = 288.938e3;
m = 200;
J = 41.35e6;
k = 6.73e5;
k_w = 1.59e6;
c = 1.43e5;
c_w = 4066;
L = 70;
S = 510.97;
Lf = L/2;
Lr = L/2; % non va cambiato?
aero_c = 8.3; %???
th = [M J m k k_w c c_w S L Lr Lf]';

d = 0;
Tend = 40;
Ts = 0.0001;
time = 0:Ts:Tend;
N = length(time);


z0 = [    0    ,    280/3.6    ,    100    ,    -25/3.6*0   ,    15*2*pi/360   ,    0      ,    300-1.9    ,    0    ,    300-1.9    ,    0    ]';
%         X          X_dot           Z         Z_dot        theta         theta_dot    z_w_r   z_w_r_dot   z_w_fr  z_w_fr_dot



z = zeros(length(z0),N);
z(:,1) = z0;
psi = zeros(1,N);
theta = zeros(1,N);
F1 = zeros(1,N);
F2 = zeros(1,N);
Vrelative_x = zeros(1,N);
Vrelative_y = zeros(1,N);
Abs_Vrelative = zeros(1,N);
u = zeros(7,N);

u(1,1) = 1e3;      % thrust force(N)
u(2:3,:) = 100;    % brake rear force(N), brake front force(N)
u(4,:) = 0.3;      % flap opening 
u(5,:) = 1;        % air-brakes opening
u(6:7,:) = 0.2;    % rear active suspension force, front active suspension force

for ind = 2:N
    z_star = z(:,ind-1) + Ts/2*aircraft(time, z(:,ind-1),u(:,ind-1),d,th);
    z(:,ind) = z(:,ind-1) + Ts/2*aircraft(time, z_star, u(:,ind-1),d,th);
    u(1,ind) = 0.9*u(1,ind-1);
    if z(3,ind) < 2 
       u(1,ind) = 0*u(1,ind-1);
    end
    [~, theta(ind), psi(ind), F1(ind), F2(ind)] = aircraft(time, z_star, u(:,ind-1),d,th);
    Abs_Vrelative(ind) = sqrt((d - z(2,ind) + aero_c*sin(z(5,ind))*z(6,ind))^2 + (-z(4,ind) - aero_c*cos(z(5,ind))*z(6,ind)^2));
    Vrelative_x(ind) = Abs_Vrelative(ind)*cos(psi(ind));
    Vrelative_y(ind) = Abs_Vrelative(ind)*sin(psi(ind));
end

figure;
plot(time, z(1,:));
xlabel('Time (s)'),ylabel('Position: X (m)')

figure;
plot(time, z(2,:));
xlabel('Time (s)'),ylabel('Horizontal velocity: X_(dot) (m)')


figure;
plot(time, z(3,:));
xlabel('Time (s)'),ylabel('Height: Z (m)')

figure;
grid on;
plot(z(1,:),z(3,:));
xlabel('Position: X (m)'),ylabel('Height: Z (m)')

figure;
plot(time, z(4,:));
xlabel('Time (s)'),ylabel('vertical velocity: Z_(dot) (m/s)')

figure;
plot(time, theta*360/2/pi);
xlabel('Time (s)'),ylabel('Pitch angle: theta (degrees)')

figure;
plot(time, psi*360/2/pi);
xlabel('Time (s)'),ylabel('Psi (degrees)')

figure;
plot(time, F1);
xlabel('Time (s)'),ylabel('Lift Force: Fl (N)')

figure;
plot(time, F2);
xlabel('Time (s)'),ylabel('Drag Force: Fd (N)')

figure;
plot(time, Abs_Vrelative);
xlabel('Time (s)'),ylabel('Absolute value of Vrel (m/s)')

figure;
plot(time, z(7,:));
xlabel('Time (s)'),ylabel('Rear wheel height (m)')

figure;
plot(time, z(9,:));
xlabel('Time (s)'),ylabel('Front wheel height (m)')

% figure;
% plot3(time, Vrelative_x, Vrelative_y);
