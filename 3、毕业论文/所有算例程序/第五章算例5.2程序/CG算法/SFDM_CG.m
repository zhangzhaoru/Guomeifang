%---------------基于CG算法的FHN模型隐式有限差分格式的求解-----------------------------
function SFDM_CG()
m1 = 256;     % x轴网格剖分节点数
m2 = 256;     % y轴网格剖分节点数
m1_half = m1/2; 
m2_half = m2/2;
N = m1_half*m2_half;
M = (m1 - 1)*(m2 - 1); % 系数矩阵A的维度

alpha1 = 2;   % x方向导数阶数
alpha2 = 2;   % y方向导数阶数
Dx = 2.5;       % x方向扩散系数
Dy = 2.5;       % y方向扩散系数
Kx = 0.0001;    % x方向扩散系数
Ky = 0.0001;    % y方向扩散系数
hx = Dx/m1;     % x方向细网格步长
hy = Dy/m2;     % y方向细网格步长
Hx = Dx/m1_half;  % x方向粗网格步长
Hy = Dy/m2_half;  % y方向粗网格步长
tau = 0.1;   % 时间步长
t = 1000;     % 终止时间
layend = t / tau;  %循环次数

%-----FHN模型中参数的选取---------%
a=0.1;        
epsilon = 0.01;
beta=0.5;
gamma = 1;
delta=0;

Ca1 = 1.0/(2 * cos(pi * alpha1 / 2)); % x方向Rieze分数阶导数定义的系数
Ca2 = 1.0/(2 * cos(pi * alpha2 / 2)); % y方向Rieze分数阶导数定义的系数

r1 = (-1) * (tau * Kx * Ca1)/(hx^alpha1);
r2 = (-1) * (tau * Ky * Ca2)/(hy^alpha2);
R1= (-1) * (tau * Kx * Ca1)/(Hx^alpha1);
R2= (-1) * (tau * Ky * Ca2)/(Hy^alpha2);

%% initial values of u,v 
u = zeros(m1-1,m2-1);
v = zeros(m1-1,m2-1);
for i = 1:m1-1
    for j = 1:m2-1
        if( ( i < m1_half ) && ( j <= m2_half ) )    % u的左下四分之一区域值为1
            u( i,j ) = 1.0;
        else
            u( i ,j  ) = 0.0;
        end
        if( ( i >= m1_half ) && ( i < m1 ) )    % v的上半区域为1
            v( i,j  ) = 0.1;
        else
            v( i,j) = 0.0;
        end
    end
end

%% initial Gruwald weights g(L)
g1 = zeros(1,m1+1);
g2 = zeros(1,m1+1);
g1(1) = 1.0;
g2(1) = 1.0;
for L = 1:m1
    g1(L+1)=-1*g1(L)*(alpha1-L+1)/L;
    g2(L+1)=-1*g2(L)*(alpha2-L+1)/L;
end

%% initial A
A=sparse(M,M);
for i = 1:m1-1
    for j = 1:m2-1
        A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j) = 1.0;  %系数矩阵对角线上的值
        if i == m1-1  % 左边界特殊处理
                A((i-1)*(m2-1)+j,(i-(1:i)+1-1)*(m2-1)+j) = A((i-1)*(m2-1)+j,(i-(1:i)+1-1)*(m2-1)+j)+(-1*r1 * g1((1:i)+1));
        else
                A((i-1)*(m2-1)+j,(i-(0:i)+1-1)*(m2-1)+j) = A((i-1)*(m2-1)+j,(i-(0:i)+1-1)*(m2-1)+j)+(-1*r1 * g1((0:i)+1));
        end
        if i == 1    % 右边界特殊处理
                A((i-1)*(m2-1)+j,(i+(1:m1-i)-1-1)*(m2-1)+j) = A((i-1)*(m2-1)+j,(i+(1:m1-i)-1-1)*(m2-1)+j)+(-1*r1 * g1((1:m1-i)+1));
        else
                A((i-1)*(m2-1)+j,(i+(0:m1-i)-1-1)*(m2-1)+j) = A((i-1)*(m2-1)+j,(i+(0:m1-i)-1-1)*(m2-1)+j)+(-1*r1 * g1((0:m1-i)+1));
        end
        if j == m2-1   % 下边界特殊处理
                A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j-(1:j)+1) = A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j-(1:j)+1)+(-1*r2* g2((1:j)+1));
        else
                A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j-(0:j)+1) = A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j-(0:j)+1)+(-1*r2 * g2((0:j)+1));
        end
        if j == 1      % 上边界特殊处理
                A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j+(1:m2-j)-1) = A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j+(1:m2-j)-1)+(-1*r2 * g2((1:m2-j)+1));
        else
                A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j+(0:m2-j)-1) = A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j+(0:m2-j)-1)+(-1*r2* g2((0:m2-j)+1));
        end
    end
        fprintf('完成矩阵A的行数：%f \n',i);     % 每完成矩阵A的一行输出一次
end
fprintf('Matrix A is over\n');      % 矩阵A计算完成
save A A;

%% 求解隐式差分格式
for lay = 1:layend
    % comute the right-hand side b
    b = u + tau * ( u .* ( 1 - u ) .* ( u - a ) - v );   % 差分格式右端项
    b=reshape(b',M,1);
    u0 = bicg(A,b);
    u=reshape(u0,m1-1,m2-1)';
    %% solve_v
    v = v + epsilon*tau*( beta * u -gamma * v -delta ); % 求解常微分方程
     
    if( mod(lay,100) ==0 )  % 每隔100层存储一次
        fprintf('时间层为：%f \n',lay);
       % eval(['save datau',num2str(lay)]);
    end    
end
contour(u)