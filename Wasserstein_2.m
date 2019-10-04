function [W_val] = Wasserstein_2(rho_0, rho_1, Nx, Nt, L, T)

dx = L / Nx;
dt = T / Nt;  % u_max*dt<=dx  
rho_0_t = rho_0/norm(rho_0,1)*Nx;
rho_1_t = rho_1/norm(rho_1,1)*Nx;

%% Solve
rhoini = kron(rho_0_t, ones(1,Nt+1));
qini = zeros(Nx+1, Nt+1);


w0 = [reshape(rhoini, Nx*(Nt+1), 1); reshape(qini, (Nx+1)*(Nt+1), 1)];
f = @(w) funJ(dx, dt, Nx, Nt, w);
options = optimoptions('fmincon','Display','iter',...
    'StepTolerance',1e-14,...
    'ConstraintTolerance',1e-14,...
    'SpecifyObjectiveGradient', true, ...
    'Checkgradients', false, 'HessianApproximation', ...
    'lbfgs', 'OptimalityTolerance', 1e-14,...
    'SpecifyConstraintGradient',true,...
    'MaxIterations', 100,...
    'MaxFunctionEvaluations', 100);
[Aeq, beq] = linfp_mat(dx, dt, Nx, Nt, rho_0_t, rho_1_t);
[A,b] = bound_mat(dx, dt, Nx, Nt);
[w fval exitflag output lambda] = fmincon(f,w0,A,b,Aeq,beq,[],[],[],options);
W_val=fval
end