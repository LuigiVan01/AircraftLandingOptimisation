function [f] = new_Flight_cost_fmincon(X,z0,nu,nz,d,Ts,Tend,ds_u,Q,R,z_ref,th)

% Parameters
N = Tend/Ts;

U  = X;

ds_T = ds_u(1,1);
ds_L = ds_u(2,1);
ds_D = ds_u(3,1);
ds_th = ds_u(4,1);

% Input sequence decomposition
fin = 0;
T = U(fin+1:fin+N/ds_T,1);
fin = fin + N/ds_T;
L = U(fin+1:fin+N/ds_L,1);
fin = fin + N/ds_L;
D = U(fin+1:fin+N/ds_D,1);
fin = fin + N/ds_D;
theta = U(fin+1:fin+N/ds_th,1);

% initialization of the cost function vector

F = zeros(nz*(N+1),1);
 
 % simulation
zsim = zeros(nz*(N+1),1);
zsim(1:nz,1) = z0;
height = zeros(N,1);
zd = zeros(nz*(N+1),1);
ztemp = z0;
u_now = [T(1,1); L(1,1); D(1,1); theta(1,1)];
cont = zeros(nu,1);
index = ones(nu,1);
 
 for ind = 2:N+1
     cont = cont + ones(nu,1); 
     % input allocation
     if cont(1,1) == ds_T+1
         index(1,1) = index(1,1) + 1;
         u_now(1,1) = T(index(1,1),1);
         cont(1,1) = 0;
     end
     
     if cont(2,1) == ds_L+1
         index(2,1) = index(2,1) + 1;
         u_now(2,1) = L(index(2,1),1);
         cont(2,1) = 0;
     end
     
     if cont(3,1) == ds_D+1
         index(3,1) = index(3,1) + 1;
         u_now(3,1) = D(index(3,1),1);
         cont(3,1) = 0;
     end
     
     if cont(4,1) == ds_th+1
         index(4,1) = index(4,1) + 1;
         u_now(4,1) = theta(index(4,1),1);
         cont(4,1) = 0;
     end
     
     % simulation
     [zdot]                              =   fly2(0,ztemp,...
                                                    u_now,d,th);
     zprime                              =   ztemp + Ts/2*zdot;
     ztemp                               =   ztemp+Ts*fly2(0,zprime,...
                                                u_now,d,th);
     zsim((ind-1)*nz+1:ind*nz,1)         =   ztemp;
     zd((ind-1)*nz+1:ind*nz,1)           =   zdot;  
     
     % Update the cost function
     F((ind-2)*nz+1:(ind-1)*nz,1)          =   Q*zdot/sqrt(N);
 end


f = F'*F;

end