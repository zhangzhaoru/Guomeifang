%---------------基于二重网格算法的FHN模型隐式有限差分格式的求解-----------------------------
function SFDM_MG()
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
        if i == m1-1  % 细网格左边界特殊处理
                A((i-1)*(m2-1)+j,(i-(1:i)+1-1)*(m2-1)+j) = A((i-1)*(m2-1)+j,(i-(1:i)+1-1)*(m2-1)+j)+(-1*r1 * g1((1:i)+1));
        else
                A((i-1)*(m2-1)+j,(i-(0:i)+1-1)*(m2-1)+j) = A((i-1)*(m2-1)+j,(i-(0:i)+1-1)*(m2-1)+j)+(-1*r1 * g1((0:i)+1));
        end
        if i == 1    % 细网格右边界特殊处理
                A((i-1)*(m2-1)+j,(i+(1:m1-i)-1-1)*(m2-1)+j) = A((i-1)*(m2-1)+j,(i+(1:m1-i)-1-1)*(m2-1)+j)+(-1*r1 * g1((1:m1-i)+1));
        else
                A((i-1)*(m2-1)+j,(i+(0:m1-i)-1-1)*(m2-1)+j) = A((i-1)*(m2-1)+j,(i+(0:m1-i)-1-1)*(m2-1)+j)+(-1*r1 * g1((0:m1-i)+1));
        end
        if j == m2-1   % 细网格下边界特殊处理
                A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j-(1:j)+1) = A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j-(1:j)+1)+(-1*r2* g2((1:j)+1));
        else
                A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j-(0:j)+1) = A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j-(0:j)+1)+(-1*r2 * g2((0:j)+1));
        end
        if j == 1      % 细网格上边界特殊处理
                A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j+(1:m2-j)-1) = A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j+(1:m2-j)-1)+(-1*r2 * g2((1:m2-j)+1));
        else
                A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j+(0:m2-j)-1) = A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j+(0:m2-j)-1)+(-1*r2* g2((0:m2-j)+1));
        end
    end
        fprintf('完成矩阵A的行数：%f \n',i);     % 每完成矩阵A的一行输出一次
end
fprintf('Matrix A is over\n');      % 矩阵A计算完成
save A A;

 %% 计算系数矩阵A1
 A1=sparse(N,N);
 for i = 1:m1_half
    for j = 1:m2_half
        A1((i-1)*m2_half+j,(i-1)*m2_half+j) = 1.0; %系数矩阵对角线上的值
        if i == m1_half    % 粗网格左边界特殊处理
                A1((i-1)*m2_half+j,(i-(1:i)+1-1)*m2_half+j) = A1((i-1)*m2_half+j,(i-(1:i)+1-1)*m2_half+j)+(-1*R1 * g1((1:i)+1));
        else
                A1((i-1)*m2_half+j,(i-(0:i)+1-1)*m2_half+j) = A1((i-1)*m2_half+j,(i-(0:i)+1-1)*m2_half+j)+(-1*R1 * g1((0:i)+1));
        end
        if i == 1          % 粗网格右边界特殊处理
                A1((i-1)*m2_half+j,(i+(1:m1_half-i)-1-1)*m2_half+j) = A1((i-1)*m2_half+j,(i+(1:m1_half-i)-1-1)*m2_half+j)+(-1*R1 * g1((1:m1_half-i)+1));
        else
                A1((i-1)*m2_half+j,(i+(0:m1_half-i)-1-1)*m2_half+j) = A1((i-1)*m2_half+j,(i+(0:m1_half-i)-1-1)*m2_half+j)+(-1*R1 * g1((0:m1_half-i)+1));
        end
        if j == m2_half    % 粗网格下边界特殊处理
                A1((i-1)*m2_half+j,(i-1)*m2_half+j-(1:j)+1) = A1((i-1)*m2_half+j,(i-1)*m2_half+j-(1:j)+1)+(-1*R2 * g2((1:j)+1));
        else
                A1((i-1)*m2_half+j,(i-1)*m2_half+j-(0:j)+1) = A1((i-1)*m2_half+j,(i-1)*m2_half+j-(0:j)+1)+(-1*R2 * g2((0:j)+1));
        end
        if j == 1          % 粗网格上边界特殊处理
                A1((i-1)*m2_half+j,(i-1)*m2_half+j+(1:m1_half-j)-1) = A1((i-1)*m2_half+j,(i-1)*m2_half+j+(1:m1_half-j)-1)+(-1*R2 * g2((1:m1_half-j)+1));
        else
                A1((i-1)*m2_half+j,(i-1)*m2_half+j+(0:m1_half-j)-1) = A1((i-1)*m2_half+j,(i-1)*m2_half+j+(0:m1_half-j)-1)+(-1*R2 * g2((0:m1_half-j)+1));
        end
    end
    fprintf('完成矩阵A1的行数：%f \n',i);     % 每完成矩阵A1的一行输出一次
 end
fprintf('Matrix A1 is over\n');      % 矩阵A1计算完成
save A1 A1;

%% 求解隐式差分格式
for lay = 1:layend
    % comute the right-hand side b
    b = u + tau * ( u .* ( 1 - u ) .* ( u - a ) - v );   % 差分格式右端项
    
   %% 松弛过程
    b=reshape(b',M,1);  %把b按行排成列向量
    u=reshape(u',M,1);  %把u按行排成列向量
    AL=tril(A,-1);      %系数矩阵A的下三角区
    AU=triu(A,1);       %系数矩阵A的上三角区
    AD=diag(diag(A));   %取A的对角线上的值构成对角阵
    for i=1:3           % Gauss_Seidel迭代3次
    u=-AD\(AL+AU)*u+AD\b;
    end
    
    %% 计算余量rh
    rh=b-A*u;
    
    %% 限制余量到粗网格上为rH
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
    %% 将校正量插值到细网格上
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
    %% 对插值后的解向量进行光滑迭代
    for i=1:3     % 雅克比迭代3次
    u=-AD\(AL+AU)*u+AD\b;
    end
    %% solve_v
    u=reshape(u,m1-1,m2-1)';
    
    if( mod(lay,100) ==0 )  % 每隔100层存储一次
        fprintf('时间层为：%f \n',lay);
        %eval(['save datau',num2str(lay)]);
    end    
    v = v + epsilon*tau*( beta * u -gamma * v -delta ); % 求解常微分方程
end
contour(u)