%%
clc;
n = 8;
E1 = 180000;
E2 = 10000;
G12 = 5000;
nu12 = 0.25;
nu21 = E2*(nu12/E1);
alpha1 = 0.2*10^-7;
alpha2 = 0.225*10^-4;
beta1 = 0;
beta2 = 0;
delta_T = 50;
delta_C = 0;
%%
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);
z = zeros(1,8);
%%
NT = zeros(3,1);
MT = zeros(3,1);
NC = zeros(3,1);
MC = zeros(3,1);
N = zeros(3,1);
M = zeros(3,1);
%%
alpha_global = zeros(3,1,n);
beta_global = zeros(3,1,n);
strain_lamina = zeros(3,1,n);
stress_lamina = zeros(3,1,n);
%%
alpha_local = [alpha1;alpha2;0];
beta_local = [beta1;beta2;0];
t = 0.125*ones(1,n);
theta = (pi/180)*[0 +45 -45 90 90 -45 +45 0];
h0 = -0.50*n*0.125;
h = zeros(1,n);
for i = 1:n
    h(i) = h0 + i*t(i);
end
%%
Q = zeros(3,3);
Q(1,1) = E1/(1 - nu12*nu21);
Q(1,2) = nu12*E2/(1 - nu12*nu21);
Q(2,2) = E2/(1 - nu12*nu21);
Q(3,3) = G12;
Q(2,1) = Q(1,2);
%%
Q_bar = zeros(3,3,n);
for i = 1:n
    Q_bar(:,:,i) = qtransform(theta(i),Q);
    if(i == 1)
        z(i) = (h(i) + h0)/2;
    else
        z(i) = (h(i) + h(i-1))/2;
    end
end
%%
for k = 1:n
    A(:,:) = A(:,:) + Q_bar(:,:,k)*t(k);
    
    if k ==1
        B(:,:) = B(:,:) + Q_bar(:,:,k)*t(k)*(h0 + h(k))/2;
        D(:,:) = D(:,:) + Q_bar(:,:,k)*(h(k)^3 - h0^3)/3;
    else
        B(:,:) = B(:,:) + Q_bar(:,:,k)*t(k)*(h(k-1) + h(k))/2;
        D(:,:) = D(:,:) + Q_bar(:,:,k)*(h(k)^3 - h(k-1)^3)/3;
    end
end
%%
for i = 1:n
    [alpha_global(:,:,i), beta_global(:,:,i)] = abtransform(alpha_local,beta_local,theta(i));
end

for i = 1:n
    NT(:,:) = NT(:,:) + delta_T*Q_bar(:,:,i)*alpha_global(:,:,i)*t(i);
    MT(:,:) = MT(:,:) + delta_T*Q_bar(:,:,i)*alpha_global(:,:,i)*t(i)*z(i);
    NC(:,:) = NC(:,:) + delta_C*Q_bar(:,:,i)*beta_global(:,:,i)*t(i);
    MC(:,:) = MC(:,:) + delta_C*Q_bar(:,:,i)*beta_global(:,:,i)*t(i)*z(i);
end
%%
Load = [(NT + NC + N);(MT + MC + M)];
Total = [A B;
         B D];

strain_global = Total\Load;
%%
for i = 1:n
    strain_lamina(:,:,i) = strain_global(1:3,:) + z(i)*strain_global(4:6,:) - alpha_global(:,:,i)*delta_T;
end

for i = 1:n
    stress_lamina(:,:,i) = qtransform(theta(i),Q_bar(:,:,i))*strain_lamina(:,:,i);
end
%%

    

