%% parameters
L = 1;
T = 1;

sigma = 0.07;


%% Grids
Nx = 300;
dx = L / Nx;
Nt = 100;
dt = T / Nt;  % u_max*dt<=dx  

%% initial/boundary conditions

%rho_0 = 2.5 * ones(Nx, 1)+2*sin(linspace(0,6*pi,Nx)');
%rho_1 = 2.5 * ones(Nx, 1)+2*cos(linspace(0,6*pi,Nx)');
r1 = 2*L/3;
r2 = L/3;
%rho_0 = linspace(1,2,Nx)'.*(2+sin(3*pi*linspace(0,1,Nx)'));
%rho_1 = linspace(1,2,Nx)'.*(2+sin(3*pi*linspace(0,1,Nx)'));

rho_0 = 1+exp(-linspace(0 - r1,L - r1,Nx).^2/2/sigma^2)/sigma;
rho_1 = 1+exp(-linspace(0 - r2,L - r2,Nx).^2/2/sigma^2)/sigma;
plot(rho_1)
