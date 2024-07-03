function [f,zsim] = Ground_cost(X,n_free,nu,d,Ts,Tend_td,Tend_b,ds_u_td,ds_u_b,Q,R,Qf,Qf_dot,x_ref,th)

N_td = Tend_td/Ts;
N_b = Tend_b/Ts;
N = N_td + N_b;
nz          = 6;

U_td = X(n_free+1:n_free+nu*N_td/ds_u_td,1);
U_b = X(n_free+nu*N_td/ds_u_td+1:n_free+nu*N_td/ds_u_td+nu*N_b/ds_u_b,1);
z0 = [0;
      X(1,1); 
      0; 
      X(2,1); 
      X(3,1);
      X(4,1)];

Nu_td=length(U_td);
Nu_b=length(U_b);

zsim        =   zeros(nz*(N+1),1);
zsim(1:nz,1) =   z0;
ztemp       =   z0;

F           =   [zeros(nz*N,1);
                zeros(nz,1);
                zeros(nz,1);
                zeros(2*(Nu_td+Nu_b)+4,1)];

u = zeros(nu,1);

for i = 1:N_td/ds_u_td
    u = U_td((i-1)*nu+1:i*nu,1);
    for ind = (i-1)*ds_u_td+2:i*ds_u_td+1
    % Update the state (RK2 simulation)
        [zdot]                              =   ground2(0,ztemp,...
                                                    u,d,th);
        zprime                              =   ztemp + Ts/2*zdot;
        ztemp                               =   ztemp+Ts*ground2(0,zprime,...
                                                    u,d,th);
        zsim((ind-1)*nz+1:ind*nz,1)           =   ztemp;

        % Update the cost function
        F((ind-2)*nz+1:(ind-1)*nz,1)          =   Q*zdot/sqrt(N);
    end
end

for j = 1:N_b/ds_u_b
    u = U_b((j-1)*nu+1:j*nu,1);
    for ind = N_td + (j-1)*ds_u_b+2:N_td + j*ds_u_b+1
    % Update the state (RK2 simulation)
        [zdot]                              =   ground2(0,ztemp,...
                                                    u,d,th);
        zprime                              =   ztemp + Ts/2*zdot;
        ztemp                               =   ztemp+Ts*ground2(0,zprime,...
                                                    u,d,th);
        zsim((ind-1)*nz+1:ind*nz,1)           =   ztemp;

        % Update the cost function
        F((ind-2)*nz+1:(ind-1)*nz,1)          =   Q*zdot/sqrt(N);
    end
end

% Update the terminal cost
F(nz*N+1:nz*(N+1),1)    =   Qf*[x_ref - zsim(N*nz+1,1); zsim(N*nz+2:nz*(N+1),1)];
F(nz*(N+1)+1:nz*(N+2),1)     =   Qf_dot*zdot;

% Update barrier function
U = [U_td; U_b];
F(nz*(N+2)+1:end,1) = barrier_input_gr(z0,U,nu);
% F(nz*(N+2)+2*Nu_td+4+1:nz*(N+2)+2*Nu_td+4+2*Nu_b+4,1) = barrier_input_gr([0;0;0;0;5*pi/180;0],U_b,nu);
f = (F'*F);
end