%---------------����CG�㷨��FHNģ����ʽ���޲�ָ�ʽ�����-----------------------------
function SFDM_CG()
m1 = 256;     % x�������ʷֽڵ���
m2 = 256;     % y�������ʷֽڵ���
m1_half = m1/2; 
m2_half = m2/2;
N = m1_half*m2_half;
M = (m1 - 1)*(m2 - 1); % ϵ������A��ά��

alpha1 = 2;   % x����������
alpha2 = 2;   % y����������
Dx = 2.5;       % x������ɢϵ��
Dy = 2.5;       % y������ɢϵ��
Kx = 0.0001;    % x������ɢϵ��
Ky = 0.0001;    % y������ɢϵ��
hx = Dx/m1;     % x����ϸ���񲽳�
hy = Dy/m2;     % y����ϸ���񲽳�
Hx = Dx/m1_half;  % x��������񲽳�
Hy = Dy/m2_half;  % y��������񲽳�
tau = 0.1;   % ʱ�䲽��
t = 1000;     % ��ֹʱ��
layend = t / tau;  %ѭ������

%-----FHNģ���в�����ѡȡ---------%
a=0.1;        
epsilon = 0.01;
beta=0.5;
gamma = 1;
delta=0;

Ca1 = 1.0/(2 * cos(pi * alpha1 / 2)); % x����Rieze�����׵��������ϵ��
Ca2 = 1.0/(2 * cos(pi * alpha2 / 2)); % y����Rieze�����׵��������ϵ��

r1 = (-1) * (tau * Kx * Ca1)/(hx^alpha1);
r2 = (-1) * (tau * Ky * Ca2)/(hy^alpha2);
R1= (-1) * (tau * Kx * Ca1)/(Hx^alpha1);
R2= (-1) * (tau * Ky * Ca2)/(Hy^alpha2);

%% initial values of u,v 
u = zeros(m1-1,m2-1);
v = zeros(m1-1,m2-1);
for i = 1:m1-1
    for j = 1:m2-1
        if( ( i < m1_half ) && ( j <= m2_half ) )    % u�������ķ�֮һ����ֵΪ1
            u( i,j ) = 1.0;
        else
            u( i ,j  ) = 0.0;
        end
        if( ( i >= m1_half ) && ( i < m1 ) )    % v���ϰ�����Ϊ1
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
        A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j) = 1.0;  %ϵ������Խ����ϵ�ֵ
        if i == m1-1  % ��߽����⴦��
                A((i-1)*(m2-1)+j,(i-(1:i)+1-1)*(m2-1)+j) = A((i-1)*(m2-1)+j,(i-(1:i)+1-1)*(m2-1)+j)+(-1*r1 * g1((1:i)+1));
        else
                A((i-1)*(m2-1)+j,(i-(0:i)+1-1)*(m2-1)+j) = A((i-1)*(m2-1)+j,(i-(0:i)+1-1)*(m2-1)+j)+(-1*r1 * g1((0:i)+1));
        end
        if i == 1    % �ұ߽����⴦��
                A((i-1)*(m2-1)+j,(i+(1:m1-i)-1-1)*(m2-1)+j) = A((i-1)*(m2-1)+j,(i+(1:m1-i)-1-1)*(m2-1)+j)+(-1*r1 * g1((1:m1-i)+1));
        else
                A((i-1)*(m2-1)+j,(i+(0:m1-i)-1-1)*(m2-1)+j) = A((i-1)*(m2-1)+j,(i+(0:m1-i)-1-1)*(m2-1)+j)+(-1*r1 * g1((0:m1-i)+1));
        end
        if j == m2-1   % �±߽����⴦��
                A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j-(1:j)+1) = A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j-(1:j)+1)+(-1*r2* g2((1:j)+1));
        else
                A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j-(0:j)+1) = A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j-(0:j)+1)+(-1*r2 * g2((0:j)+1));
        end
        if j == 1      % �ϱ߽����⴦��
                A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j+(1:m2-j)-1) = A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j+(1:m2-j)-1)+(-1*r2 * g2((1:m2-j)+1));
        else
                A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j+(0:m2-j)-1) = A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j+(0:m2-j)-1)+(-1*r2* g2((0:m2-j)+1));
        end
    end
        fprintf('��ɾ���A��������%f \n',i);     % ÿ��ɾ���A��һ�����һ��
end
fprintf('Matrix A is over\n');      % ����A�������
save A A;

%% �����ʽ��ָ�ʽ
for lay = 1:layend
    % comute the right-hand side b
    b = u + tau * ( u .* ( 1 - u ) .* ( u - a ) - v );   % ��ָ�ʽ�Ҷ���
    b=reshape(b',M,1);
    u0 = bicg(A,b);
    u=reshape(u0,m1-1,m2-1)';
    %% solve_v
    v = v + epsilon*tau*( beta * u -gamma * v -delta ); % ��ⳣ΢�ַ���
     
    if( mod(lay,100) ==0 )  % ÿ��100��洢һ��
        fprintf('ʱ���Ϊ��%f \n',lay);
       % eval(['save datau',num2str(lay)]);
    end    
end
contour(u)