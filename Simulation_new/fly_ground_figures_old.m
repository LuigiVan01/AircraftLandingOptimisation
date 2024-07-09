clear
close all
clc

%% Parameters
M = 22e3;
m = 130;
J = 100e3;
k_f = 6.73e5;
k_r = 1.59e4;            % modified: different (smaller) value for the front stiffness
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
Tend = 170;
Ts = 0.01;
time = 0:Ts:Tend;
z_iniz = 500;

V_cruise = sqrt((2*M*g)/(rho*S*cl0)) + d;

N = length(time);
z0 = [0, V_cruise, z_iniz, 0, 0, 0, z_iniz-2, 0, z_iniz-2, 0]';
z = zeros(length(z0),N);
z(:,1) = z0;
psi = zeros(1,N);
theta = zeros(1,N);
psi_a = zeros(1,N);
alpha_out = zeros(1,N);
F1 = zeros(1,N);   % Fs_r
F2 = zeros(1,N);   % Fs_fr
F3 = zeros(1,N);   % Fl
F4 = zeros(1,N);   % Fd
F5 = zeros(1,N);
check = zeros(1,N);
Vrelative_x = zeros(1,N);
Vrelative_y = zeros(1,N);
Abs_Vrelative = zeros(1,N);
u = zeros(8,N);
u(1,1) = 0.6 %Thrust
u(2:3,:) = 8000; %Freni
u(4,:) = 0; %uL
u(5,:) = 0; %uD
u(6:7,:) = 1e3; %Sospensioni
u(8,:) = 0;
g = 9.81;
flag = 0;

% tic
% for ind = 2:N
%     if z(3,ind-1) > 1 && flag == 0
%         z_star = z(:,ind-1) + Ts/2*fly(time, z(:,ind-1),u(:,ind-1),d,th);
%         z(:,ind) = z(:,ind-1) + Ts/2*fly(time, z_star, u(:,ind-1),d,th);
%         u(1,ind) = 0.02*u(1,1);
%         [~, theta(ind), psi(ind), psi_a(ind), alpha_out(ind), F1(ind), F2(ind), F3(ind), F4(ind), F5(ind)] = ...
%             fly(time, z_star, u(:,ind-1),d,th);        
%         Abs_Vrelative(ind) = sqrt((d - z(2,ind) + aero_c*sin(z(5,ind))*z(6,ind))^2 + ... 
%             (-z(4,ind) - aero_c*cos(z(5,ind))*z(6,ind))^2);
%     %     Abs_Vrelative(ind) = sqrt((d - z(2,ind))^2 + (-z(4,ind))^2);
%         Vrelative_x(ind) = Abs_Vrelative(ind)*cos(psi(ind));
%         Vrelative_y(ind) = Abs_Vrelative(ind)*sin(psi(ind));
%     else
%         flag = 1;
%         z_star = z(:,ind-1) + Ts/2*ground(time, z(:,ind-1),u(:,ind-1),d,th);
%         z(:,ind) = z(:,ind-1) + Ts/2*ground(time, z_star, u(:,ind-1),d,th);
%         u(1,ind) = 0.02*u(1,1);
%         [~, theta(ind), psi(ind), psi_a(ind), alpha_out(ind), F1(ind), F2(ind), F3(ind), F4(ind), F5(ind)] = ...
%             ground(time, z_star, u(:,ind-1),d,th);        
%         Abs_Vrelative(ind) = sqrt((d - z(2,ind) + aero_c*sin(z(5,ind))*z(6,ind))^2 + ... 
%             (-z(4,ind) - aero_c*cos(z(5,ind))*z(6,ind))^2);
%     %     Abs_Vrelative(ind) = sqrt((d - z(2,ind))^2 + (-z(4,ind))^2);
%         Vrelative_x(ind) = Abs_Vrelative(ind)*cos(psi(ind));
%         Vrelative_y(ind) = Abs_Vrelative(ind)*sin(psi(ind));
%     end
% end
% toc

%Simulation with ODE45
tic
for ind = 2:N
    if z(3,ind-1) > 0 && flag == 0
        z_out_temp = ode45(@(t,z)fly2(t, z, u(:,ind-1), d, th), [0 Ts], z(:,ind-1));
        z(:,ind) = z_out_temp.y(:,end);
%         u(1,ind) = 0.02*u(1,ind-1);
        [~, theta(ind), psi(ind), psi_a(ind), alpha_out(ind), F1(ind), F2(ind), F3(ind), F4(ind), F5(ind)] = ...
        fly2(time, z(:,ind-1), u(:,ind-1),d,th);        
%         Abs_Vrelative(ind) = sqrt((d - z(2,ind) + aero_c*sin(z(5,ind))*z(6,ind))^2 + ... 
%         (-z(4,ind) - aero_c*cos(z(5,ind))*z(6,ind))^2);
        Abs_Vrelative(ind) = sqrt((d - z(2,ind))^2 + (-z(4,ind))^2);
        Vrelative_x(ind) = Abs_Vrelative(ind)*cos(psi(ind));
        Vrelative_y(ind) = Abs_Vrelative(ind)*sin(psi(ind));
        
        % input management flight zone
        if ind*Ts > 10
            u(4,ind) = 1;  %lift 
            u(5,ind) = 0;  %drag
            u(2:3,ind) = 0; %brakes
            u(6:7,ind) = 0; %suspensions
            if z(3,ind-1) > 14.3
                u(1,ind) = u(1,ind-1) - 10;
                check(ind) = F5(ind)/tan(-3*pi/180) - V_cruise*M/(ind*Ts - 10) - F3(ind)*cos(177*pi/180 - pi/2) ...
                    - F4(ind)*cos(177*pi/180);
                if u(1,ind) <  check(ind)
                    u(1,ind) = check(ind);
                end
                if u(1,ind) < 3000 
                    u(1,ind) = 3000;
                end
                u(8,ind) = 0*pi/180; %pitch angle
            else 
%                 u(1,ind) = 0.95*u(1,ind-100); %thrust
                u(8,ind) = 10*pi/180; %pitch angle
            end
        else
            u(1,ind) = -(F3(ind)*cos(psi(ind) - pi/2) + F4(ind)*cos(psi(ind)))/cos(theta(ind));
        end
    else
        flag = 1;
        z_out_temp = ode45(@(t,z)ground2(t, z, u(:,ind-1), d, th), [0 Ts], z(:,ind-1));
        z(:,ind) = z_out_temp.y(:,end);
        u(1,ind) = 0.02*u(1,ind-1);
        [~, theta(ind), psi(ind), psi_a(ind), alpha_out(ind), F1(ind), F2(ind), F3(ind), F4(ind), F5(ind)] = ...
        ground2(time, z(:,ind-1), u(:,ind-1),d,th);        
        Abs_Vrelative(ind) = sqrt((d - z(2,ind) + aero_c*sin(z(5,ind))*z(6,ind))^2 + ... 
        (-z(4,ind) - aero_c*cos(z(5,ind))*z(6,ind))^2);
    %   Abs_Vrelative(ind) = sqrt((d - z(2,ind))^2 + (-z(4,ind))^2);
        Vrelative_x(ind) = Abs_Vrelative(ind)*cos(psi(ind));
        Vrelative_y(ind) = Abs_Vrelative(ind)*sin(psi(ind));
        
        %input management ground zone
        u(4,ind) = 0;  %lift 
        u(5,ind) = 0;  %drag
        if u(1,ind-1) > 0   %thrust
            u(1,ind) = u(1,ind-1) - 50;
        else
            u(1, ind) = 0;
        end
        %pitch not used
        % 2. Active suspension
        if z(4,ind) < 0
            u(6:7,ind) = 0;
        elseif z(4,ind) > 0
            u(6:7,ind) = u(6:7,ind-1) + 1e3;
        else
            u(6:7, ind) = 0;
        end
        %brakes
        if z(2,ind) < 0.5
            u(2:3, ind) = 0;
        end
            
    end

end
toc

% Simulation with forward euler
% tic
% for ind = 2:N
%     z(:,ind) = z(:,ind-1) + Ts*aircraft_correction(time, z(:,ind-1), u(:,ind-1),d, th);
%     u(1,ind) = 0.02*u(1,ind-1);
% end
% toc
%% Plots
figure(1); % x,y position wrt time
subplot 211
plot(time, z(1,:));
xlabel('Time [s]','Interpreter','latex');
ylabel('X [m]','Interpreter','latex');
title('X,Z positions wrt time','Interpreter','latex');
grid;
subplot 212
plot(time, z(3,:));
grid;
xlabel('Time [s]','Interpreter','latex');
ylabel('Z [m]','Interpreter','latex');


figure(2); %trajectory
plot(z(1,:),z(3,:));
grid;
xlabel('X [m]','Interpreter','latex');
ylabel('Z [m]','Interpreter','latex');
title('Trajectory','Interpreter','latex');

figure(3);
subplot 221
plot(time, theta*360/2/pi,'b');
grid;
xlabel('Time [s]','Interpreter','latex');
ylabel('Deg','Interpreter','latex')
title('Pitch Angle $\vartheta$','Interpreter','latex');

subplot 222
plot(time, psi*360/2/pi,'k');
grid;
xlabel('Time [s]','Interpreter','latex');
ylabel('Deg','Interpreter','latex')
title('$\Psi$','Interpreter','latex');

subplot 223
plot(time, psi_a*360/2/pi,'r');
grid;
xlabel('Time [s]','Interpreter','latex');
ylabel('Deg','Interpreter','latex')
title('$\Psi_{a}$','Interpreter','latex');

subplot 224
plot(time, alpha_out*360/2/pi);
grid;
xlabel('Time [s]','Interpreter','latex');
ylabel('Deg','Interpreter','latex')
title('Attack Angle $\alpha$','Interpreter','latex');

figure(4); % x,y,vrel speed profiles
plot(time, z(4,:),'b');
hold on;
grid;
plot(time, z(2,:)),'k';
plot(time, Abs_Vrelative,'r');
xlabel('Time [s]','interpreter','latex');
ylabel('Velocity [$ms^{-1}$]','Interpreter','latex');
title('Speed profiles','Interpreter','latex');
legend('$\dot{Z}$','$\dot{X}$','$V_{rel}$','interpreter','latex');

figure(5); %drag and lift forces
plot(time, F3, 'r');
hold on;
grid;
plot(time, F4, 'b');
xlabel('Time [s]','Interpreter','latex');
ylabel('Force [N]','Interpreter','latex');
title('Drag and Lift Forces','Interpreter','latex');
legend('$F_{Lift}$','$F_{Drag}$','interpreter','latex');

figure(6) %lift analysis
subplot 211
plot(time, psi_a,'r');
xlabel('Time [s]','Interpreter','latex');
ylabel('[deg]','Interpreter','latex');
title('Drag and Lift parameters analysis','interpreter','latex');
hold on;
grid;
plot(time, psi);
legend('$\Psi_{a}$','$\Psi$','interpreter','latex');

subplot 212
plot(time, Abs_Vrelative,'b');
grid;
xlabel('Time [s]','Interpreter','latex');
ylabel('Velocity [$ms^{-1}$]','Interpreter','latex');

figure(7)
title('ul, ud, T Inputs','Interpreter','latex');
grid;
hold on
plot(time, u(4,:));
plot(time, u(5,:));
plot(time, u(1,:)/10000);
legend('$u_{L}$','$u_{D}$','$T$','interpreter','latex');
xlabel('Time [s]','interpreter','latex')
hold off

figure(8)
subplot 211
plot(time, z(5,:)*180/pi);
grid;
xlabel('Time [s]','Interpreter','Latex');
ylabel('[deg]','Interpreter','latex');
title('Pitch Angle $\vartheta$','Interpreter','latex');

subplot 212
plot(time, z(6,:)*180/pi);
grid;
xlabel('Time [s]','Interpreter','latex');
ylabel('[$degs^{-1}$]','Interpreter','latex');
title('Pitching Velocity $\dot{\vartheta}$','Interpreter','latex');

figure(9); 
title('Elastic Forces');
hold on
plot(time, F1, 'r');
plot(time, F2), 'b';
legend({'y = Fl','y = Fd'});
hold off

figure(10)
title('Total force acting on z axis');
plot(time, F5/M);   % over M to compute the acceleration


% figure(11)
% plot(time, check); 
