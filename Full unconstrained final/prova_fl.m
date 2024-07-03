close all
z = zeros(6,N+1);
zd = zeros(6,N+1);
z(:,1) = z0;
time = 1:N+1;
nu = 4;
nz = 6;
u = zeros(4,N/ds_u);
u(1,:) = 0.6;
u(2,:) = 0.5;
u(3,:) = 0;
u(4,:) = 0;
ztemp = z0;

Xstar = zeros(4*N/ds_u,1);
for j = 1:N/ds_u
    for i=1:4
        Xstar((j-1)*4+i,1) = u(i,j);
    end
end

% RK2 simulation
        
for i = 1:N/ds_u
    for ind = (i-1)*ds_u+2:i*ds_u+1
%         if z(3,ind-1) < 25 
%             Xstar((i-1)*4+4,1) = 3/20;
%             Xstar((i-1)*4+2,1) = 0.4;
% %             Xstar((i-1)*4+1,1) = 0.3;
%         end
        if z(3,ind-1) < 15
            Xstar((i-1)*4+4,1) = 6/20;
            Xstar((i-1)*4+2,1) = 0.5;
%             Xstar((i-1)*4+1,1) = 0.2;
        end
        if z(3,ind-1) < 1
            Xstar((i-1)*4+4,1) = 0/20;
            Xstar((i-1)*4+2,1) = 1;
            Xstar((i-1)*4+3,1) = 0.8;
%             Xstar((i-1)*4+1,1) = 0.3;
        end
        u_check = Xstar((i-1)*nu+1:i*nu,1);
        u(:,i) = u_check;
        [zdot]                              =   fly2(0,ztemp,u_check,d,th);
        zd(:,ind) = zdot;
        zprime                              =   ztemp + Ts/2*zdot;
        ztemp                               =   ztemp+Ts*fly2(0,zprime,u_check,d,th);
        z(:,ind)           =   ztemp;
    end
end


F = Flight_cost_GN(Xstar,z0,d,Ts,Tend,ds_u,Q,R,Qf,z_ref,th);
F_state = F(1:nz*N,1);
F_final_state = F(nz*N+1:nz*(N+1),1);
F_barrier = F(nz*(N+2)+1:end,1);
[~,z_sim] = Flight_cost_GN(Xstar,z0,d,Ts,Tend,ds_u,Q,R,Qf,z_ref,th);
zsim = zeros(nz,N+1);
for ind = 1:N+1
    zsim(:,ind) = z_sim((ind-1)*nz+1:ind*nz,1);
end

for j=1:6
    figure(j)
    hold on
    if j == 5 || j == 6
        plot(time,180/pi*z(j,:));
        plot(time,180/pi*zsim(j,:), 'b*');
        plot(time,180/pi*zd(j,:));
    else
        plot(time,z(j,:));
        plot(time,zsim(j,:), 'b*');
        plot(time,zd(j,:));
    end
    hold off
end
for j=1:4
    figure(6+j)
    if j == 4
        plot(20*u(j,:));
    else
        plot(u(j,:));
    end
end
figure(11)
hold on
plot(z(1,:), z(3,:));
plot(zsim(1,:), zsim(3,:),'b*');
hold off

% save initialguess_fl.mat Xstar
