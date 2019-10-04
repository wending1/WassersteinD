function mov = plotfig(L, Nx, T, Nt, rho, u)
    mov(1:Nt+1) = struct('cdata', [], 'colormap', []);
    dx = L / Nx;
    dt = T / Nt;
    for n = 1 : Nt
        plot(linspace(dx/2,L-dx/2,Nx), rho(:,n),'b',...
            linspace(0,L-dx,Nx), u(:,n),'r');
        title("time = "+num2str((n-1)*dt,'%-1.2f'));
        legend('density', 'velocity');
        ylim([min(min(u))-0.2,max(max(max(u)),max(max(rho)))+0.2]);
        pause(0.1)
        mov(n) = getframe();
end