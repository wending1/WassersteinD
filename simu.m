%% parameters

L = 1;
T = 1;

sigma = 0.07;


%% Grids
Nx = 300;
dx = L / Nx;
Nt = 300;
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

rho_0 = rho_0/norm(rho_0,1)*Nx;
rho_1 = rho_1/norm(rho_1,1)*Nx;

%% Solve
rhoini = kron(rho_0, ones(1,Nt+1));
qini = zeros(Nx+1, Nt+1);
%rho0 = 0.2*ones(Nx, Nt+1);
%q0 = 0.5 * ones(Nx, Nt+1);
%rho0 = 1.5+rand(Nx, Nt+1);
%q0 = rand(Nx, Nt+1);

w0 = [reshape(rhoini, Nx*(Nt+1), 1); reshape(qini, (Nx+1)*(Nt+1), 1)];
f = @(w) funJ(dx, dt, Nx, Nt, w);
options = optimoptions('fmincon','Display','iter',...
    'StepTolerance',1e-15,...
    'ConstraintTolerance',1e-15,...
    'SpecifyObjectiveGradient', true, ...
    'Checkgradients', false, 'HessianApproximation', ...
    'lbfgs', 'OptimalityTolerance', 1e-15,...
    'SpecifyConstraintGradient',true,...
    'MaxIterations', 100,...
    'MaxFunctionEvaluations', 100);
[Aeq, beq] = linfp_mat(dx, dt, Nx, Nt, rho_0, rho_1);
[A,b] = bound_mat(dx, dt, Nx, Nt);
[w fval exitflag output lambda] = fmincon(f,w0,A ,b,Aeq,beq,[],[],[],options);