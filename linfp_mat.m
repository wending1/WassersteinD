function [Aeq, beq] = linfp_mat(dx, dt, Nx, Nt, rho_0, rho_1)
    %Aeq = sparse(Nx*(Nt+1), Nx*(Nt+1)+(Nx+1)*(Nt+1));
    beq = zeros(Nx*Nt+2*Nx+2*(Nt+1), 1);
    
    index = 0;
    eqs = 0;
    idx = zeros(4*Nt*Nx+2*Nx+2*(Nt+1),1);
    idy = zeros(4*Nt*Nx+2*Nx+2*(Nt+1),1);
    idv = zeros(4*Nt*Nx+2*Nx+2*(Nt+1),1);
    
    for j = 1 : Nx
        eqs = eqs+1;
        index = index + 1;
        idx(index) = eqs;
        idy(index) = labelrho(Nx, j, 1);
        idv(index) = 1;
        beq(eqs) = rho_0(j);
        
        eqs = eqs+1;
        index = index + 1;
        idx(index) = eqs;
        idy(index) = labelrho(Nx, j, Nt+1);
        idv(index) = 1;
        beq(eqs) = rho_1(j);
    end
    for n = 1 : Nt+1
        eqs = eqs+1;
        index = index + 1;
        idx(index) = eqs;
        idy(index) = labelq(Nx, Nt, 1, n);
        idv(index) = 1;
        
        eqs = eqs+1;
        index = index + 1;
        idx(index) = eqs;
        idy(index) = labelq(Nx, Nt, Nx+1, n);
        idv(index) = 1;
    end
    
    
    for n = 1 : Nt
    for j = 1 : Nx
        eqs = eqs +1;
        
        index = index + 1;
        idx(index) = eqs;
        idy(index) = labelrho(Nx, j, n+1);
        idv(index) = 1;
        
        index = index + 1;
        idx(index) = eqs;
        idy(index) = labelrho(Nx, j, n);
        idv(index) = -1;
        
        index = index + 1;
        idx(index) = eqs;
        idy(index) = labelq(Nx, Nt, j+1, n);
        idv(index) = dt/dx;
        
        index = index + 1;
        idx(index) = eqs;
        idy(index) = labelq(Nx, Nt, j, n);
        idv(index) = -dt/dx;
    end
    end
    
    Aeq = sparse(idx, idy, idv, Nt*Nx+2*Nx+2*(Nt+1), Nx*(Nt+1)+(Nx+1)*(Nt+1));
end

function l = labelrho(Nx, j, n)
    l = (n-1) * Nx + j;
end

function l = labelq(Nx, Nt, j, n)
    l = Nx * (Nt+1);
    l = l + (n-1) * (Nx + 1) + j;
end