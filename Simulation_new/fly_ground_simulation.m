clear all
close all
clc

%% Model parameters

M = 22e3;
m = 130;
J = 100e3;
k_f = 6.73e5;
k_r = 1.59e4;           % modified: different (smaller) value for the front stiffness
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


g   =  9.81;
d   =  0;
nz  =  6;


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

T0 = 0.5;
L0 = 1;
D0 = 1;
th0 = 0.1;

u_fl = [T0; L0; D0; th0];

%% Inputs ground
T0 = 0;
L0 = 0.1;
D0 = 1;
B0 = 1;
Far0 = 0;
Faf0 =0;

u_gr = [T0; L0; D0; B0; Far0; Faf0];


%% RK2 Simulation

zsim = zeros(nz,N+1);
zsim(:,1) = z0;
zd = zeros(nz,N+1);
ztemp = z0;
flag = 0;
 
tic
 for ind = 2:N+1
     if ztemp(3) > 1 && flag == 0

        [zdot]              =   fly2(0,ztemp,u_fl,d,th);

        zprime              =   ztemp + Ts/2*zdot;
        ztemp               =   ztemp + Ts*fly2(0,zprime,u_fl,d,th);

        zsim(:,ind)         =   ztemp;
        zd(:,ind)           =   zdot;

        %u_fl(1) = 0.02*u_fl(1);

     else
        flag=1;
     
        [zdot]              =   ground2(0,ztemp,u_gr,d,th);

        zprime              =   ztemp + Ts/2*zdot;
        ztemp               =   ztemp + Ts*ground2(0,zprime,u_gr,d,th);

        zsim(:,ind)         =   ztemp;
        zd(:,ind)           =   zdot;

     end
 end
t_RK2=toc;

%% ode45 Simulation 

zsim_ode45 = zeros(nz,N+1);
zsim_ode45(:,1) = z0;
zd_ode45 = zeros(nz,N+1);
ztemp_ode45 = z0;
flag = 0;

%u_fl(1) = 0.5;
tic
for ind=2:N+1
    if zsim_ode45(3,ind-1) > 1 && flag == 0

        ztemp_ode45                =   ode45(@(t,z)fly2(t,z,u_fl,...
                                       0,th),[0 Ts], zsim_ode45(:,ind-1));

        zsim_ode45(:,ind)          =   ztemp_ode45.y(:,end);


       [zdot_ode45]                =   fly2(0,ztemp_ode45.y(:,end),u_fl,d,...
                                       th);

        zd_ode45(:,ind)            =   zdot_ode45;

        %u_fl(1) = 0.02*u_fl(1);
    else

        flag=1;

        ztemp_ode45                 =   ode45(@(t,z)ground2(t,z,u_gr,...
                                        0,th),[0 Ts],zsim_ode45(:,ind-1));

        zsim_ode45(:,ind)           =   ztemp_ode45.y(:,end);

        [zdot_ode45]                =   ground2(0,ztemp_ode45.y(:,end),...
                                        u_gr,d,th);

        zd_ode45(:,ind)             =   zdot_ode45;
    end
end
t_o45=toc;


%% Plot of the States

time = 0:Ts:Tend_fl;

figure(1)
hold on
plot(time,zsim(1,:),'b','DisplayName','position RK2');
plot(time,zd(1,:),'r','DisplayName','speed RK2');
plot(time,zsim_ode45(1,:),'g','DisplayName','position ode45');
plot(time,zd_ode45(1,:),'c','DisplayName','speed ode45');
grid
legend
title('Horizontal position overview',"Interpreter","Latex")
xlabel('Time',"Interpreter","Latex");
hold off

figure(2)
hold on
plot(time,zsim(2,:),'b','DisplayName','speed RK2');
plot(time,zd(2,:),'r','DisplayName','acceleration RK2');
plot(time,zsim_ode45(2,:),'g','DisplayName','speed ode45');
plot(time,zd_ode45(2,:),'c','DisplayName','acceleretion ode45');
grid 
legend
title('Horizontal speed overview',"Interpreter","Latex")
xlabel('Time',"Interpreter","Latex");
hold off

figure(3)
hold on
plot(time,zsim(3,:),'b','DisplayName','position RK2');
plot(time,zd(3,:),'r','DisplayName','speed RK2');
plot(time,zsim_ode45(3,:),'g','DisplayName','position ode45');
plot(time,zd_ode45(3,:),'c','DisplayName','speed ode45');
grid 
legend
title('Vertical position overview',"Interpreter","Latex")
xlabel('Time',"Interpreter","Latex");
hold off

figure(4)
hold on
plot(time,zsim(4,:),'b','DisplayName','speed RK2');
plot(time,zd(4,:),'r','DisplayName','acceleration RK2');
plot(time,zsim_ode45(4,:),'g','DisplayName','speed ode45');
plot(time,zd_ode45(4,:),'c','DisplayName','acceleretion ode45');

grid
legend
title('Vertical acceleration overview',"Interpreter","Latex")
xlabel('Time',"Interpreter","Latex");
hold off

figure(5)
hold on
plot(time,180/pi*zsim(5,:),'b','DisplayName','position RK2');
plot(time,180/pi*zd(5,:),'r','DisplayName','speed RK2');
plot(time,180/pi*zsim_ode45(5,:),'g','DisplayName','position ode45');
plot(time,180/pi*zd_ode45(5,:),'c','DisplayName','speed ode45');
grid 
legend
title('Pitch overview',"Interpreter","Latex")
xlabel('Time',"Interpreter","Latex");
hold off

figure(6)
hold on
plot(time,180/pi*zsim(6,:),'b','DisplayName','speed RK2');
plot(time,180/pi*zd(6,:),'r','DisplayName','acceleration RK2');
plot(time,180/pi*zsim_ode45(6,:),'g','DisplayName','speed ode45');
plot(time,180/pi*zd_ode45(6,:),'c','DisplayName','acceleretion ode45');
grid 
legend
title('Pitch acceleration overview',"Interpreter","Latex")
xlabel('Time',"Interpreter","Latex");
hold off

figure(7)
hold on
plot(zsim(1,:), zsim(3,:),'b','DisplayName','Trajectory RK2');
plot(zsim_ode45(1,:), zsim_ode45(3,:),'r','DisplayName','Trajectory ode45');
legend
grid
title('Trajectory',"Interpreter","Latex")