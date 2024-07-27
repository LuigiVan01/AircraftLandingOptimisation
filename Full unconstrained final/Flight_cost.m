function [f,zsim] = Flight_cost(X,z0,d,Ts,Tend,ds_u,Q,R,Qf,z_ref,th) %#codegen

U = X;
nz = length(z0);
nu = 4;
Nu=length(U);
N = Tend/Ts;

zsim        =   zeros(nz*(N+1),1);
height = zeros(N-1,1);
zsim(1:nz,1) =   z0;
ztemp       =   z0;

F           =   [zeros(nz*N,1);
                zeros(nz,1);
                zeros(2*Nu+N-1,1)];

for i = 1:N/ds_u
    u = U((i-1)*nu+1:i*nu,1);
    for ind = (i-1)*ds_u+2:i*ds_u+1
    % Update the state (RK2 simulation)
        [zdot]                              =   fly2(0,ztemp,...
                                                    u,d,th);
        zprime                              =   ztemp + Ts/2*zdot;
        ztemp                               =   ztemp+Ts*fly2(0,zprime,...
                                                    u,d,th);
        zsim((ind-1)*nz+1:ind*nz,1)          =   ztemp;
        if ind < N+1
            height(ind-1) = ztemp(3,1);
        end
        % Update the cost function
        F((ind-2)*nz+1:(ind-1)*nz,1)          =   Q*zdot/sqrt(N);
    end
end

% Update the terminal cost
F(nz*N+1:nz*N+nz,1)    =   Qf*(zsim(N*nz+1:nz*(N+1),1) - z_ref(:,1));

%Update barrier function
F(nz*(N+1)+1:end,1) = barrier_input_fl(U,nu,height);

f = (F'*F);