function SFDM_CG(N)
% N 为空间网格点数
m1 = N;
m2 = N;
M = (m1 - 1)*(m2 - 1);
alpha1 =1.5;
alpha2 = 1.5;
Dx = 1;
Dy = 1;
Kx = 1;
Ky = 1;
t = 1;
b = zeros(m1-1,m2-1);
u = zeros(m1-1,m2-1);
g1 = zeros(1,m1+1);
g2 = zeros(1,m1+1);
hx = Dx/m1; 
hy = Dy/m2;
tau = 1/100;

% Phi1t=1-exp(-tau);
% Phi2t=(1-exp(-epsilon*tau))/epsilon;
Phi1t=tau;

Ca1 = -1.0/(2 * cos(pi * alpha1 / 2));
Ca2 = -1.0/(2 * cos(pi * alpha2 / 2));

r1 = (Phi1t * Kx * Ca1)/(hx^alpha1);
r2 = (Phi1t * Ky * Ca2)/(hy^alpha2);

layend = t / Phi1t;
x=linspace(0,Dx,m1-1);
y=linspace(0,Dy,m2-1);
x1=repmat(x,1,m1-1)';
x=repmat( x',1,m1-1);
y=repmat(y,m2-1,1);
y1=reshape(y,M,1);
%% initial u
u=10*x.^2.*(1-x).^2.*y.^2.*(1-y).^2;
uexact=10*exp(-t)*x.^2.*(1-x).^2.*y.^2.*(1-y).^2;
% figure(1)
% mesh(uexact);
%% initial g(L)
g1(1) = 1.0;
g2(1) = 1.0;
for L = 1:m1
    g1(L+1)=-1*g1(L)*(alpha1-L+1)/L;
    g2(L+1)=-1*g2(L)*(alpha2-L+1)/L;
end
 %load('A.mat')
%  load('A1.mat')
%% initial A
A=sparse(M,M);
for i = 1:m1-1
    for j = 1:m2-1
        A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j) = 1.0;
        if i == m1-1
                A((i-1)*(m2-1)+j,(i-(1:i)+1-1)*(m2-1)+j) = A((i-1)*(m2-1)+j,(i-(1:i)+1-1)*(m2-1)+j)+(-1*r1 * g1((1:i)+1));
        else
                A((i-1)*(m2-1)+j,(i-(0:i)+1-1)*(m2-1)+j) = A((i-1)*(m2-1)+j,(i-(0:i)+1-1)*(m2-1)+j)+(-1*r1 * g1((0:i)+1));
        end
        if i == 1
                A((i-1)*(m2-1)+j,(i+(1:m1-i)-1-1)*(m2-1)+j) = A((i-1)*(m2-1)+j,(i+(1:m1-i)-1-1)*(m2-1)+j)+(-1*r1 * g1((1:m1-i)+1));
        else
                A((i-1)*(m2-1)+j,(i+(0:m1-i)-1-1)*(m2-1)+j) = A((i-1)*(m2-1)+j,(i+(0:m1-i)-1-1)*(m2-1)+j)+(-1*r1 * g1((0:m1-i)+1));
        end
        if j == m2-1
                A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j-(1:j)+1) = A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j-(1:j)+1)+(-1*r2* g2((1:j)+1));
        else
                A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j-(0:j)+1) = A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j-(0:j)+1)+(-1*r2 * g2((0:j)+1));
        end
        if j == 1
                A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j+(1:m2-j)-1) = A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j+(1:m2-j)-1)+(-1*r2 * g2((1:m2-j)+1));
        else
                A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j+(0:m2-j)-1) = A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j+(0:m2-j)-1)+(-1*r2* g2((0:m2-j)+1));
        end
    end
        fprintf('完成矩阵A的行数：%f \n',i);
end
fprintf('Matrix A is over\n');
save A A;

for lay = 1:layend
   % comute the right-hand side b
   Fu=-u-u.*u;
   f1=f(x,y,lay*tau,alpha1,alpha2,Kx,Ky,Ca1,Ca2);
   b=u+Phi1t*Fu+Phi1t*f1;  
    
    b=reshape(b',M,1);
    u0 = bicg(A,b);
    u=reshape(u0,m1-1,m2-1)';
    % 输出 techplot 格式数据文件
    if( mod(lay,200) ==0 )
        fprintf('时间层为：%f \n',lay);
    end
end
% figure(2)
% mesh(u)
% figure(3)
format long
u_error_temp=uexact-u;
[u_error_temp_n,u_error_temp_m] = size(u_error_temp);
u1_temp = u_error_temp(1:u_error_temp_n-1,1:u_error_temp_m-1);
u2_temp = u_error_temp(2:u_error_temp_n  ,1:u_error_temp_m-1);
u3_temp = u_error_temp(1:u_error_temp_n-1,2:u_error_temp_m  );
u4_temp = u_error_temp(2:u_error_temp_n  ,2:u_error_temp_m  );
u_error = (u1_temp+u2_temp+u3_temp+u4_temp)/4;
error = sqrt(sum(sum(u_error.^2*hx*hy)));
error2 =norm(uexact-u,2);
error_inf =max(max(abs(uexact-u)));
end
% mesh(uexact-u)
% u1