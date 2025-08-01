function [b]  = barrier_input_gr(z0,u,nu)

theta_in = z0(5,1);
theta_dot_in = z0(6,1);
Nu = length(u); 
u_barr = zeros(nu,Nu);

theta_max = 10/180*pi;
theta_min = 2/180*pi;
theta_dot_max = 10/180*pi;
theta_dot_min = -5/180*pi;
T_max = 1;
T_min = 0;
ul_max = 1;
ul_min = 0;
ud_max = 1;
ud_min = 0;
Brake_max = 1;
Brake_min = 0;
Fa_r_max = 1;
Fa_r_min = -1;
Fa_fr_max = 1;
Fa_fr_min = -1;

alpha_T = 0.001;
alpha_ul = 0.001;
alpha_ud = 0.001;
alpha_brake = 0.001;
alpha_Fa_r = 0.001;
alpha_Fa_fr = 0.001;
alpha_theta = 0.00001;
beta_T = 7e2;
beta_ul = 7e2;
beta_ud = 7e2;
beta_brake = 7e2;
beta_Fa_r = 7e2;
beta_Fa_fr = 7e2;
beta_theta = 11e2;

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
   
   barrier_cost_max(4,ind) = alpha_brake*exp(-beta_brake*(Brake_max - u_barr(4,ind)));
   barrier_cost_min(4,ind) = alpha_brake*exp(-beta_brake*(u_barr(4,ind) - Brake_min));
   
   barrier_cost_max(5,ind) = alpha_Fa_r*exp(-beta_Fa_r*(Fa_r_max - u_barr(5,ind)));
   barrier_cost_min(5,ind) = alpha_Fa_r*exp(-beta_Fa_r*(u_barr(5,ind) - Fa_r_min));
   
   barrier_cost_max(6,ind) = alpha_Fa_fr*exp(-beta_Fa_fr*(Fa_fr_max - u_barr(6,ind)));
   barrier_cost_min(6,ind) = alpha_Fa_fr*exp(-beta_Fa_fr*(u_barr(6,ind) - Fa_fr_min));
end

b = zeros(2*Nu+4,1);
for h = 1:Nu/nu
    for j=1:nu
        b((h-1)*nu+j,1) = barrier_cost_max(j,h);
        b(Nu+(h-1)*nu+j,1) = barrier_cost_min(j,h);
    end
end
b(2*Nu+1) = alpha_theta*exp(-beta_theta*(theta_max - theta_in));
b(2*Nu+2) = alpha_theta*exp(-beta_theta*(theta_in - theta_min));
b(2*Nu+3) = alpha_theta*exp(-beta_theta*(theta_dot_max - theta_dot_in));
b(2*Nu+4) = alpha_theta*exp(-beta_theta*(theta_dot_in - theta_dot_min));

end