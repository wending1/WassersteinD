function [A, b] = bound_mat(dx, dt, Nx, Nt)
    A = sparse(Nx*(Nt+1), Nx*(Nt+1)+(Nx+1)*(Nt+1));
    b = zeros(Nx*(Nt+1), 1);
    cnt = 0;
    for n = 1 : Nt+1
    for j = 1 : Nx
        cnt = cnt + 1;
        A(cnt, labelrho(Nx, j, n)) = -1;
    end
    end
%     for n = 1 : Nt+1
%     for j = 1 : Nx
%         cnt = cnt + 1;
%         A(cnt, labelrho(Nx, j, n)) = 1;
%         b(cnt) = rho_max;
%     end
%     end
%     for n = 1 : Nt+1
%     for j = 1 : Nx
%         cnt = cnt + 1;
%         A(cnt, labelq(Nx, Nt, j, n)) = -1;
%         b(cnt) = 0;
%     end
%     end
%     for n = 1 : Nt+1
%     for j = 1 : Nx
%         cnt = cnt + 1;
%         A(cnt, labelq(Nx, Nt, j+1, n)) = 1;
%         A(cnt, labelrho(Nx, j, n)) = -u_max;
%         b(cnt) = 0;
%     end
%     end
end

function l = labelrho(Nx, j, n)
    l = (n-1) * Nx + j;
end