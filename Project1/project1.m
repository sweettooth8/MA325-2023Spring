clc;
clear;
h = 2^(-4);
N = 1/h - 1;
x = (0:N+1)*h;
y = x;
u = @(x,y) x^2*y^2*(1-x)^2*(1-y)^2;% exact solution
f = @(x,y) (12*x^2-12*x+2)*y^2*(1-y)^2+x^2*(1-x)^2*(12*y^2-12*y+2)-100*x^2*y^2*(1-x)^2*(1-y)^2;
%%  Use the five-point central finite difference method to discretize above
D_1 = zeros(N,N);
for i = 1:N
    D_1(i,i) = -2;
end
for j = 1:N-1
    D_1(j,j+1) = 1;
    D_1(j+1,j) = 1;
end
D_1 = D_1./h^2;
A = kron(eye(N),D_1) + kron(D_1,eye(N)) - 100 * eye(N*N);
F = zeros(N*N,1);
for i = 1:N
    for j = 1:N
        F(j+(i-1)*N,1) = f(x(1,j+1),y(1,i+1));
    end
end

%% Solve the linear system by Discrete Sine Transform.
lam = eig(A);
sol_dst = dst2(idst2(F)./lam); % sol_dst is the solution of the linear system by Discrete Sine Transform
%% SOR iteration method.
w = 2 / (1+sin(pi/N)); %optimal relaxation parameter
x_0 = zeros(N*N,1); % initial guess
sol_sor = sor(A,F,x_0,w,100); % sol_sor is the solution of the linear system by SOR iteration method
%% Compute the error in L2 norm.
U = zeros(N*N,1); % exact solution in vector
for i = 1:N
    for j = 1:N
        U(j+(i-1)*N,1) = u(x(1,j+1),y(1,i+1));
    end
end
E_dst = U - sol_dst;
E_sor = U - sol_sor;
norm_2_dst = norm(E_dst); % the error of dst in L2 norm
norm_2_sor = norm(E_sor); % the error of sor in L2 norm
%% Compute the error in L-Inf norm
norm_inf_dst = norm(E_dst,"inf"); % the error of dst in L-inf norm
norm_inf_sor = norm(E_sor,"inf"); % the error of sor in L-inf norm
%% plot
X = [2^(-4) 2^(-5) 2^(-6) 2^(-7)];
% dst in L2
Y_dst_2 = [0.023816185742856 0.050003918771661 0.101196584643166 0.202980689678211];
% dst in L-inf
Y_dst_inf = [0.003638923914462  0.0038406351980220 0.101196584643166 0.073460260454717];
% sor in L2
Y_sor_L2 = [1.116544021167664e-04 0.001752577889570 0.014681254589841 0.003902268226593];
% sor in L-inf
Y_sor_inf = [2.466006771741970e-05 1.294496013043694e-04 5.352148701918069e-04 0.001402458643383];

subplot(2,2,1)
loglog(X,Y_dst_2)
hold on
loglog(X,10000*X.^2)
legend("error","ch^2")
title("the dst error in L2 norm")
subplot(2,2,2)
loglog(X,Y_dst_inf)
hold on
loglog(X,10000*X.^2)
legend("error","ch^2")
title("the dst error in L-max norm")
subplot(2,2,3)
loglog(X,Y_sor_L2)
hold on
loglog(X,10000*X.^2)
legend("error","ch^2")
title("the sor error in L2 norm")
subplot(2,2,4)
loglog(X,Y_sor_inf)
hold on
loglog(X,10000*X.^2)
legend("error","ch^2")
title("the sor error in L-max norm")
