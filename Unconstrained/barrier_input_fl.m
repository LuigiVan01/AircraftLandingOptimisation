function [b]  = barrier_input_fl(u,nu,height)
Nu = length(u);
N = length(height);
u_barr = zeros(nu,Nu);
T_max = 1;
T_min = 0;
ul_max = 1;
ul_min = 0;
ud_max = 1;
ud_min = 0;
theta_max = 1;
theta_min = -1;

alpha_T = 0.01;
alpha_ul = 0.01;
alpha_ud = 0.01;
alpha_theta = 0.01;
beta_T = 8e2;
beta_ul = 8e2;
beta_ud = 8e2;
beta_theta = 8e2;

alpha = 0.00001;
beta = 3e2;

% umax - u > 0;
% u - umin > 0;

Nu = length(u);
for ind = 1:Nu/nu
    u_barr(:,ind) = u((ind-1)*nu+1:ind*nu,1);    
end

barrier_cost_max = zeros(nu,Nu);
barrier_cost_min = zeros(nu,Nu);

for ind = 1:Nu/nu
   barrier_cost_max(1,ind) = alpha_T*exp(-beta_T*(T_max - u_barr(1,ind)));
   barrier_cost_min(1,ind) = alpha_T*exp(-beta_T*(u_barr(1,ind) - T_min));
   
   barrier_cost_max(2,ind) = alpha_ul*exp(-beta_ul*(ul_max - u_barr(2,ind)));
   barrier_cost_min(2,ind) = alpha_ul*exp(-beta_ul*(u_barr(2,ind) - ul_min));
   
   barrier_cost_max(3,ind) = alpha_ud*exp(-beta_ud*(ud_max - u_barr(3,ind)));
   barrier_cost_min(3,ind) = alpha_ud*exp(-beta_ud*(u_barr(3,ind) - ud_min));
   
   barrier_cost_max(4,ind) = alpha_theta*exp(-beta_theta*(theta_max - u_barr(4,ind)));
   barrier_cost_min(4,ind) = alpha_theta*exp(-beta_theta*(u_barr(4,ind) - theta_min));
end

b = zeros(2*Nu+N,1);
for h = 1:Nu/nu
    for j=1:nu
        b((h-1)*nu+j,1) = barrier_cost_max(j,h)/sqrt(Nu);
        b(Nu+(h-1)*nu+j,1) = barrier_cost_min(j,h)/sqrt(Nu);
    end
end
for i=1:N
    b(2*Nu+i,1) = alpha*exp(-beta*(height(i) - 0))/sqrt(N);
end

end