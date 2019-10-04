function [val, grad] = funJ(dx, dt, Nx, Nt, w)
    rho = reshape(w(1:Nx*(Nt+1)), Nx, Nt+1);
    q = reshape(w(Nx*(Nt+1)+1:Nx*(Nt+1)+(Nx+1)*(Nt+1)), Nx+1, Nt+1);
    val = J(dx, dt, Nx, Nt, rho, q);
    if nargout > 1
        [gradrho, gradq] = gradJ(dx, dt, Nx, Nt, rho, q);
        grad = [reshape(gradrho, Nx*(Nt+1), 1);...
                reshape(gradq, (Nx+1)*(Nt+1), 1)];
        
    end
end
    