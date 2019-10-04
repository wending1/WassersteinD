function [gradrho, gradq] = gradJ(dx, dt, Nx, Nt, rho, q)
    gradrho = zeros(Nx, Nt+1);
    gradrho(:, 2:Nt) = - (q(1:Nx,2:Nt).^2 + q(2:Nx+1,2:Nt).^2) / 4 ./ rho(:,2:Nt).^2;
    gradrho(:, 1) = -(q(1:Nx,1).^2 + q(2:Nx+1,1).^2) / 8 ./ rho(:,1).^2;
    gradrho(:, Nt+1) = -(q(1:Nx,Nt+1).^2 + q(2:Nx+1,Nt+1).^2) / 8 ./ rho(:,Nt+1).^2;
    gradrho = gradrho * dx * dt;
    
    gradq = zeros(Nx+1, Nt+1);
    gradq(2:Nx, 2:Nt) = 0.5 * q(2:Nx, 2:Nt) .* ...
        (1 ./ rho(1:Nx-1,2:Nt) + 1 ./ rho(2:Nx,2:Nt));
    gradq(1, 2:Nt) = 0.5 * q(1, 2:Nt) ./ rho(1,2:Nt);
    gradq(Nx+1, 2:Nt) = 0.5 * q(Nx+1, 2:Nt) ./ rho(Nx,2:Nt);
    
    gradq(2:Nx, 1) = 0.5 / 2 * q(2:Nx, 1) .* ...
        (1 ./ rho(1:Nx-1,1) + 1 ./ rho(2:Nx,1));
    gradq(1, 1) = 0.5 /2  * q(1, 1) ./ rho(1, 1);
    gradq(Nx+1, 1) = 0.5 / 2 * q(Nx+1, 1) ./ rho(Nx, 1);
    
    gradq(2:Nx, Nt+1) = 0.5 / 2 * q(2:Nx, Nt+1) .* ...
        (1 ./ rho(1:Nx-1,Nt+1) + 1 ./ rho(2:Nx, Nt+1));
    gradq(1, Nt+1) = 0.5 /2  * q(1, Nt+1) ./ rho(1, Nt+1);
    gradq(Nx+1, Nt+1) = 0.5 / 2 * q(Nx+1, Nt+1) ./ rho(Nx, Nt+1);
    gradq = gradq * dt * dx;
end