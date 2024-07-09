clear all
close all
clc

%% Model parameters

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
nz =   6;


%% Parameters
Ts          =   0.01;                        % Sampling time
Tend_fl     =   15;  
N = Tend_fl/Ts;


%% Initial condition for the Flight operation

hor_speed_in = 100;
height_in = 50;
vert_speed_in = -5;
pitch_in = 0;
pitch_dot_in = 0;

z0 = [0; hor_speed_in; height_in; vert_speed_in; pitch_in; pitch_dot_in];

%% Inputs fligh

T0 = 0.6;
L0 = 1;
D0 = 1;
th0 = 0.1;

u_fl = [T0; L0; D0; th0];

%% Inputs ground
T0 = 0;
L0 = 0.1;
D0 = 0.5;
B0 = 1;
Far0 = 0;
Faf0 =0;

u_gr = [T0; L0; D0; B0; Far0; Faf0];


%% Simulation

zsim = zeros(nz*(N+1),1);
zsim(1:nz,1) = z0;
height = zeros(N,1);
zd = zeros(nz*(N+1),1);
ztemp = z0;

 
 for ind = 2:N+1
     [zdot]                              =   fly2(0,ztemp,...
                                                    u_fl,d,th);
     zprime                              =   ztemp + Ts/2*zdot;
     ztemp                               =   ztemp+Ts*fly2(0,zprime,...
                                                u_fl,d,th);
     zsim((ind-1)*nz+1:ind*nz,1)         =   ztemp;
     zd((ind-1)*nz+1:ind*nz,1)           =   zdot;
     if ind > 9*N/10 && ind ~= N+1      % the reason for that is explained in the flight main script
         height(ind-1,1)                       =   ztemp(3,1);
     end

     % if ztemp(3)<1
     %     for i=ind:N+1
     %           [zdot]                              =   ground2(0,ztemp,...
     %                                                u_gr,d,th);
     %            zprime                              =   ztemp + Ts/2*zdot;
     %            ztemp                               =   ztemp+Ts*ground2(0,zprime,...
     %                                            u_gr,d,th);
     %            zsim((i-1)*nz+1:i*nz,1)         =   ztemp;
     %            zd((i-1)*nz+1:i*nz,1)           = zdot;
     %     end
     %     ind=N+1;
     % end
     % 
 end

N_fl=N;

z = zeros(nz, N_fl+1);
z_d = zeros(nz, N_fl+1);

for ind = 1:N_fl+1
    z(:,ind) = zsim((ind-1)*nz+1:ind*nz);
    z_d(:,ind) = zd((ind-1)*nz+1:ind*nz);
end



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