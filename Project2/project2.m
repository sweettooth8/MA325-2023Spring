

clc;
clear;
N = 255;
h = 2*pi/(N+1);
dt = 0.01;
T = 10;
tn = T/dt;
e = 0.1;
b = 2/e;
x = (0:(N+1))*h;
y = x';
uin = 0.05*(2*rand(N,N)-1);
aver = sum(sum(uin))/(N^2);
uin = uin - aver;
U = zeros(N,N,T/dt);
U(:,:,1) = uin;
[xx,yy] = meshgrid(x(2:end-1),y(2:end-1));
lam = (2*(cos(xx)+cos(yy)-2))./h^2;

x_i = linspace(0,2*pi,N+2);
x_i = x_i(2:end-1);
y_i = x_i;
[X,Y] = meshgrid(x_i,y_i);
energy = zeros(1,T/dt);

for i = 1:T/dt
    F_i = (1+dt*b+dt/e).*U(:,:,i) - (dt/e)*(U(:,:,i).^3);
    A = (1+dt*b)*ones(N)-dt*e*lam;
    U(:,:,i+1) = real(ifft2(fft2(F_i)./A));
    [dx,dy] = gradient(U(:,:,i),x_i,y_i);
    f = (1/(4*e))*(U(:,:,i).^2-1).^2+(e/2)*(dx.^2+dy.^2);
    energy(1,i)=trapz(y_i,trapz(x_i,f,2));


end
subplot(2,3,1);
mesh(xx,yy,uin)
axis off;
title("T = 0")

subplot(2,3,2)
mesh(xx,yy,U(:,:,41))
axis off
title("T = 0.4")

subplot(2,3,3)
mesh(xx,yy,U(:,:,101))
axis off
title("T = 1")

subplot(2,3,4)
mesh(xx,yy,U(:,:,401))
axis off
title("T = 4")

subplot(2,3,5)
mesh(xx,yy,U(:,:,1001))
axis off
title("T = 10")

subplot(2,3,6)
plot(dt:dt:T,energy)
xlabel("时间T")
ylabel("能量E(u)")
title("the energy evolution against time")
