L = 1;
T = 1;

sigma = 0.07;
r1 = 2*L/3;
r2 = L/3;
clear wvalue;

%% Grids
Nx = 300;
dx = L / Nx;
Nt = 300;
dt = T / Nt;  % u_max*dt<=dx  
i=1;
rho_1 = 1+exp(-linspace(0 - r2,L - r2,Nx).^2/2/sigma^2)/sigma;
rho_1 = rho_1/norm(rho_1,1)*Nx;
for r=L/5:0.02:L*2/5


rho_0 = 1+exp(-linspace(0 - r,L - r,Nx).^2/2/sigma^2)/sigma;
rho_0 = rho_0+rho_1;

rho_0 = rho_0/norm(rho_0,1)*Nx;


wvalue(i)=Wasserstein_2(rho_0, rho_1, Nx, Nt, L, T);
i=i+1;
end


