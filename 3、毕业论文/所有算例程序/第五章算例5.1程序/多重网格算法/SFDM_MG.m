function SFDM_MG(N)
% N 为空间网格节点数
m1 = N;
m2 = N;
m1_half = m1/2-1;
m2_half = m2/2-1;
N=m1_half*m2_half;
M = (m1 - 1)*(m2 - 1);
alpha1 = 1.5;
alpha2 = 1.5;
Dx = 1;
Dy = 1;
tau = 0.001;
Kx = 1;
Ky = 1;
t = 1;
b = zeros(m1-1,m2-1);
u = zeros(m1-1,m2-1);
g1 = zeros(1,m1+1);
g2 = zeros(1,m1+1);
hx = Dx/m1; %细网格步长
hy = Dy/m2; 
Hx = Dx/m1_half; %粗网格步长
Hy = Dy/m2_half;
Phi1t=tau;

Ca1 = 1.0/(2 * cos(pi * alpha1 / 2));
Ca2 = 1.0/(2 * cos(pi * alpha2 / 2));

r1 = (-1) * (Phi1t* Kx * Ca1)/(hx^alpha1);
r2 = (-1) * (Phi1t * Ky * Ca2)/(hy^alpha2);

rr1= (-1) * (Phi1t* Kx * Ca1)/(Hx^alpha1);
rr2= (-1) * (Phi1t * Ky * Ca2)/(Hy^alpha2);

layend = t / tau;
x=linspace(0,Dx,m1-1);
y=linspace(0,Dy,m2-1);
x1=repmat(x,1,m1-1)';
x=repmat(x',1,m1-1);
y=repmat(y,m2-1,1);
y1=reshape(y,M,1);

%% initial u
u=10*x.^2.*(1-x).^2.*y.^2.*(1-y).^2;
uexact=10*exp(-t)*x.^2.*(1-x).^2.*y.^2.*(1-y).^2;
figure(1)
mesh(uexact);
%% initial g(L)
g1(1) = 1.0;
g2(1) = 1.0;
gg=0;
for L = 1:m1
    g1(L+1)=-1*g1(L)*(alpha1-L+1)/L;
    g2(L+1)=-1*g2(L)*(alpha2-L+1)/L;
    gg=gg+g1(L);
end
%   load('A.mat')
%   load('A1.mat')
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

 %% 计算系数矩阵A1
 A1=sparse(N,N);
 for i = 1:m1_half
    for j = 1:m2_half
        A1((i-1)*m2_half+j,(i-1)*m2_half+j) = 1.0;
        if i == m1_half
                A1((i-1)*m2_half+j,(i-(1:i)+1-1)*m2_half+j) = A1((i-1)*m2_half+j,(i-(1:i)+1-1)*m2_half+j)+(-1*rr1 * g1((1:i)+1));
        else
                A1((i-1)*m2_half+j,(i-(0:i)+1-1)*m2_half+j) = A1((i-1)*m2_half+j,(i-(0:i)+1-1)*m2_half+j)+(-1*rr1 * g1((0:i)+1));
        end
        if i == 1
                A1((i-1)*m2_half+j,(i+(1:m1_half-i)-1-1)*m2_half+j) = A1((i-1)*m2_half+j,(i+(1:m1_half-i)-1-1)*m2_half+j)+(-1*rr1 * g1((1:m1_half-i)+1));
        else
                A1((i-1)*m2_half+j,(i+(0:m1_half-i)-1-1)*m2_half+j) = A1((i-1)*m2_half+j,(i+(0:m1_half-i)-1-1)*m2_half+j)+(-1*rr1 * g1((0:m1_half-i)+1));
        end
        if j == m2_half
                A1((i-1)*m2_half+j,(i-1)*m2_half+j-(1:j)+1) = A1((i-1)*m2_half+j,(i-1)*m2_half+j-(1:j)+1)+(-1*rr2 * g2((1:j)+1));
        else
                A1((i-1)*m2_half+j,(i-1)*m2_half+j-(0:j)+1) = A1((i-1)*m2_half+j,(i-1)*m2_half+j-(0:j)+1)+(-1*rr2 * g2((0:j)+1));
        end
        if j == 1
                A1((i-1)*m2_half+j,(i-1)*m2_half+j+(1:m1_half-j)-1) = A1((i-1)*m2_half+j,(i-1)*m2_half+j+(1:m1_half-j)-1)+(-1*rr2 * g2((1:m1_half-j)+1));
        else
                A1((i-1)*m2_half+j,(i-1)*m2_half+j+(0:m1_half-j)-1) = A1((i-1)*m2_half+j,(i-1)*m2_half+j+(0:m1_half-j)-1)+(-1*rr2 * g2((0:m1_half-j)+1));
        end
    end
    fprintf('完成矩阵A1的行数：%f \n',i);
end
save A1 A1;
fprintf('Matrix A1 is over\n');

for lay = 1:layend
   % comute the right-hand side b
    Fu=-u-u.*u;
    f1=f(x,y,lay*tau,alpha1,alpha2,Kx,Ky,Ca1,Ca2);
    b=u+Phi1t*Fu+Phi1t*f1;   
    b=reshape(b',M,1);
    
    u=reshape(u',M,1); %把u按行排成列向量
    AL=tril(A,-1);%系数矩阵A的下三角区
    AU=triu(A,1);%系数矩阵A的上三角区
    AD=diag(diag(A));
    for i=1:5
    u=-AD\(AL+AU)*u+AD\b;
    end
    p=[1/4 1/2 1/4;1/2 1 1/2;1/4 1/2 1/4];
    R=1/2*p';
    %% 计算余量rh
    rh=b-A*u;
    %% 限制余量到二重网格上为rH
    rhh=reshape(rh,m1-1,m2-1)';
    rH=zeros(m1_half,m2_half);
    for i=1:m1_half
        for j=1:m2_half
            if i==1&&j==1 %左下角点
                rH(i,j)=4/9*(rhh(i,j)+1/2*(rhh(i,j+1)+rhh(i+1,j))+1/4*rhh(i+1,j+1));
            elseif i==1&&j==m2_half  %右下角点
                rH(i,j)=4/9*(rhh(i,2*j-1)+1/2*(rhh(i,2*j-2)+rhh(i+1,2*j-1))+1/4*rhh(i+1,2*j-2));
            elseif i==m1_half&&j==1  %左上角点
                rH(i,j)=4/9*(rhh(2*i-1,j)+1/2*(rhh(2*i-2,j)+rhh(2*i-1,j+1))+1/4*rhh(2*i-2,j+1));
            elseif i==m1_half&&j==m2_half %右上角点
                 rH(i,j)=4/9*(rhh(2*i-1,2*j-1)+1/2*(rhh(2*i-1,2*j-2)+rhh(2*i-2,2*j-1))+1/4*rhh(2*i-2,2*j-2));
            elseif i==1 %下边界
                rH(i,j)=1/3*(rhh(i,2*j-1)+1/2*(rhh(i,2*j-2)+rhh(i,2*j)+rhh(i+1,2*j-1))+1/4*(rhh(i+1,2*j-2)+rhh(i+1,2*j)));
            elseif j==1 %左边界
                rH(i,j)=1/3*(rhh(2*i-1,j)+1/2*(rhh(2*i-1,j+1)+rhh(2*i-2,j)+rhh(2*i,j))+1/4*(rhh(2*i-2,j+1)+rhh(2*i,j+1)));
            elseif i==m1_half %上边界
                rH(i,j)=1/3*(rhh(2*i-1,2*j-1)+1/2*(rhh(2*i-2,2*j-1)+rhh(2*i-1,2*j-2)+rhh(2*i-1,2*j))+1/4*(rhh(2*i-2,2*j-2)+rhh(2*i-2,2*j)));
            elseif j==m2_half %右边界
                rH(i,j)=1/3*(rhh(2*i-1,2*j-1)+1/2*(rhh(2*i-2,2*j-1)+rhh(2*i-1,2*j-2)+rhh(2*i,2*j-1))+1/4*(rhh(2*i-2,2*j-2)+rhh(2*i,2*j-2)));
            else  %内点
                rH(i,j)=1/4*(rhh(2*i-1,2*j-1)+1/2*(rhh(2*i-2,2*j-1)+rhh(2*i-1,2*j-2)+rhh(2*i,2*j-1)+rhh(2*i-1,2*j))+1/4*(rhh(2*i-2,2*j-2)+rhh(2*i-2,2*j)+rhh(2*i,2*j-2)+rhh(2*i,2*j)));
            end
        end
    end
    %% 在粗网格上精确计算A1*eH=rH
    rHH=reshape(rH',N,1); %把矩阵转换为向量
    eH=A1\rHH;
    eHH=reshape(eH,m1_half,m2_half)';
    eh=zeros(m1-1,m2-1);
    %将校正量插值到细网格上
    for i=1:m1_half
        for j=1:m2_half
            if i<m1_half&&j<m2_half
                eh(2*i-1,2*j-1)=eHH(i,j);
                eh(2*i-1,2*j)=1/2*(eHH(i,j)+eHH(i,j+1));
                eh(2*i,2*j-1)=1/2*(eHH(i,j)+eHH(i+1,j));  
                eh(2*i,2*j)=1/4*(eHH(i,j)+eHH(i,j+1)+eHH(i+1,j)+eHH(i+1,j+1));
            elseif i==m1_half&&j<m2_half
                    eh(2*i-1,2*j-1)=eHH(i,j);
                    eh(2*i-1,2*j)=1/2*(eHH(i,j)+eHH(i,j+1));
            elseif j==m2_half&&i<m1_half
                    eh(2*i-1,2*j-1)=eHH(i,j);
                    eh(2*i,2*j-1)=1/2*(eHH(i,j)+eHH(i+1,j)); 
            end
        end
    end
    ehh=reshape(eh',M,1);
    u=u+ehh;
    %对插值后的解向量进行光滑迭代，雅克比迭代5次
    for i=1:5
    u=-AD\(AL+AU)*u+AD\b;
    end
    %%solve_v
    u=reshape(u,m1-1,m2-1)';
    if( mod(lay,100) ==0 )
        fprintf('时间层为：%f \n',lay);
    end
end
figure(2)
mesh(u)
figure(3)
u_error_temp=uexact-u;
[u_error_temp_n,u_error_temp_m] = size(u_error_temp);
u1_temp = u_error_temp(1:u_error_temp_n-1,1:u_error_temp_m-1);
u2_temp = u_error_temp(2:u_error_temp_n  ,1:u_error_temp_m-1);
u3_temp = u_error_temp(1:u_error_temp_n-1,2:u_error_temp_m  );
u4_temp = u_error_temp(2:u_error_temp_n  ,2:u_error_temp_m  );
u_error = (u1_temp+u2_temp+u3_temp+u4_temp)/4;
error = sqrt(sum(sum(u_error.^2*hx*hy)))
error2 =norm(uexact-u,2)
error_inf =norm(uexact-u,inf)
mesh(uexact-u)
% u1