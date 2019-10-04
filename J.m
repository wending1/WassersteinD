function val = J(dx, dt, Nx, Nt, rho, q)
    %M = zeros(Nx, Nt+1);
    M = 0.5 * ((q(1:Nx,:).^2 + q(2:Nx+1,:).^2) / 2 ./ rho(1:Nx,:));

    val = dx * dt * sum(sum(M(:,2:Nt)));
    val = val + dx * dt / 2 * (sum(M(:, 1)) + sum(M(:, Nt+1)));
end