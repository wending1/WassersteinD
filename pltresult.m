 rho = reshape(w(1:Nx*(Nt+1)), Nx, Nt+1);
 q = reshape(w(Nx*(Nt+1)+1:Nx*(Nt+1)+(Nx+1)*(Nt+1)), Nx+1, Nt+1);

%{
for i =1 :Nt+1
    
    
    x=linspace(0,1,Nx);
    for j = 1:Nx
        y(j) = rho(j, i);
    end
    plot(x,y);
    hold on
end
 %}
x=linspace(0,L,Nx);
y = linspace(0,L, Nx);
for n = 1 : Nt+1
    for j = 1:Nx
        y(j) = rho(j, n);
    end
    plot(x,y);
        title("time = "+num2str((n-1)*dt,'%-1.2f'));
        legend('rho');

        %pause(0.001)
        mov(n) = getframe();
end
