%---------------���ڶ��������㷨��FHNģ����ʽ���޲�ָ�ʽ�����-----------------------------
function FHN_MG()
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
        if i == m1-1  % ϸ������߽����⴦��
                A((i-1)*(m2-1)+j,(i-(1:i)+1-1)*(m2-1)+j) = A((i-1)*(m2-1)+j,(i-(1:i)+1-1)*(m2-1)+j)+(-1*r1 * g1((1:i)+1));
        else
                A((i-1)*(m2-1)+j,(i-(0:i)+1-1)*(m2-1)+j) = A((i-1)*(m2-1)+j,(i-(0:i)+1-1)*(m2-1)+j)+(-1*r1 * g1((0:i)+1));
        end
        if i == 1    % ϸ�����ұ߽����⴦��
                A((i-1)*(m2-1)+j,(i+(1:m1-i)-1-1)*(m2-1)+j) = A((i-1)*(m2-1)+j,(i+(1:m1-i)-1-1)*(m2-1)+j)+(-1*r1 * g1((1:m1-i)+1));
        else
                A((i-1)*(m2-1)+j,(i+(0:m1-i)-1-1)*(m2-1)+j) = A((i-1)*(m2-1)+j,(i+(0:m1-i)-1-1)*(m2-1)+j)+(-1*r1 * g1((0:m1-i)+1));
        end
        if j == m2-1   % ϸ�����±߽����⴦��
                A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j-(1:j)+1) = A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j-(1:j)+1)+(-1*r2* g2((1:j)+1));
        else
                A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j-(0:j)+1) = A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j-(0:j)+1)+(-1*r2 * g2((0:j)+1));
        end
        if j == 1      % ϸ�����ϱ߽����⴦��
                A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j+(1:m2-j)-1) = A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j+(1:m2-j)-1)+(-1*r2 * g2((1:m2-j)+1));
        else
                A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j+(0:m2-j)-1) = A((i-1)*(m2-1)+j,(i-1)*(m2-1)+j+(0:m2-j)-1)+(-1*r2* g2((0:m2-j)+1));
        end
    end
        fprintf('��ɾ���A��������%f \n',i);     % ÿ��ɾ���A��һ�����һ��
end
fprintf('Matrix A is over\n');      % ����A�������
save A A;

 %% ����ϵ������A1
 A1=sparse(N,N);
 for i = 1:m1_half
    for j = 1:m2_half
        A1((i-1)*m2_half+j,(i-1)*m2_half+j) = 1.0; %ϵ������Խ����ϵ�ֵ
        if i == m1_half    % ��������߽����⴦��
                A1((i-1)*m2_half+j,(i-(1:i)+1-1)*m2_half+j) = A1((i-1)*m2_half+j,(i-(1:i)+1-1)*m2_half+j)+(-1*R1 * g1((1:i)+1));
        else
                A1((i-1)*m2_half+j,(i-(0:i)+1-1)*m2_half+j) = A1((i-1)*m2_half+j,(i-(0:i)+1-1)*m2_half+j)+(-1*R1 * g1((0:i)+1));
        end
        if i == 1          % �������ұ߽����⴦��
                A1((i-1)*m2_half+j,(i+(1:m1_half-i)-1-1)*m2_half+j) = A1((i-1)*m2_half+j,(i+(1:m1_half-i)-1-1)*m2_half+j)+(-1*R1 * g1((1:m1_half-i)+1));
        else
                A1((i-1)*m2_half+j,(i+(0:m1_half-i)-1-1)*m2_half+j) = A1((i-1)*m2_half+j,(i+(0:m1_half-i)-1-1)*m2_half+j)+(-1*R1 * g1((0:m1_half-i)+1));
        end
        if j == m2_half    % �������±߽����⴦��
                A1((i-1)*m2_half+j,(i-1)*m2_half+j-(1:j)+1) = A1((i-1)*m2_half+j,(i-1)*m2_half+j-(1:j)+1)+(-1*R2 * g2((1:j)+1));
        else
                A1((i-1)*m2_half+j,(i-1)*m2_half+j-(0:j)+1) = A1((i-1)*m2_half+j,(i-1)*m2_half+j-(0:j)+1)+(-1*R2 * g2((0:j)+1));
        end
        if j == 1          % �������ϱ߽����⴦��
                A1((i-1)*m2_half+j,(i-1)*m2_half+j+(1:m1_half-j)-1) = A1((i-1)*m2_half+j,(i-1)*m2_half+j+(1:m1_half-j)-1)+(-1*R2 * g2((1:m1_half-j)+1));
        else
                A1((i-1)*m2_half+j,(i-1)*m2_half+j+(0:m1_half-j)-1) = A1((i-1)*m2_half+j,(i-1)*m2_half+j+(0:m1_half-j)-1)+(-1*R2 * g2((0:m1_half-j)+1));
        end
    end
    fprintf('��ɾ���A1��������%f \n',i);     % ÿ��ɾ���A��һ�����һ��
 end
fprintf('Matrix A1 is over\n');      % ����A�������
save A1 A1;

%% �����ʽ��ָ�ʽ
for lay = 1:layend
    % comute the right-hand side b
    b = u + tau * ( u .* ( 1 - u ) .* ( u - a ) - v );   % ��ָ�ʽ�Ҷ���
    
   %% �ɳڹ���
    b=reshape(b',M,1);  %��b�����ų�������
    u=reshape(u',M,1);  %��u�����ų�������
    AL=tril(A,-1);      %ϵ������A����������
    AU=triu(A,1);       %ϵ������A����������
    AD=diag(diag(A));   %ȡA�ĶԽ����ϵ�ֵ���ɶԽ���
    for i=1:3           % Gauss_Seidel����3��
    u=-AD\(AL+AU)*u+AD\b;
    end
    
    %% ��������rh
    rh=b-A*u;
    
    %% ������������������ΪrH
    rhh=reshape(rh,m1-1,m2-1)';
    rH=zeros(m1_half,m2_half);
    for i=1:m1_half
        for j=1:m2_half
            if i==1&&j==1 %���½ǵ�
                rH(i,j)=4/9*(rhh(i,j)+1/2*(rhh(i,j+1)+rhh(i+1,j))+1/4*rhh(i+1,j+1));
            elseif i==1&&j==m2_half  %���½ǵ�
                rH(i,j)=4/9*(rhh(i,2*j-1)+1/2*(rhh(i,2*j-2)+rhh(i+1,2*j-1))+1/4*rhh(i+1,2*j-2));
            elseif i==m1_half&&j==1  %���Ͻǵ�
                rH(i,j)=4/9*(rhh(2*i-1,j)+1/2*(rhh(2*i-2,j)+rhh(2*i-1,j+1))+1/4*rhh(2*i-2,j+1));
            elseif i==m1_half&&j==m2_half %���Ͻǵ�
                 rH(i,j)=4/9*(rhh(2*i-1,2*j-1)+1/2*(rhh(2*i-1,2*j-2)+rhh(2*i-2,2*j-1))+1/4*rhh(2*i-2,2*j-2));
            elseif i==1 %�±߽�
                rH(i,j)=1/3*(rhh(i,2*j-1)+1/2*(rhh(i,2*j-2)+rhh(i,2*j)+rhh(i+1,2*j-1))+1/4*(rhh(i+1,2*j-2)+rhh(i+1,2*j)));
            elseif j==1 %��߽�
                rH(i,j)=1/3*(rhh(2*i-1,j)+1/2*(rhh(2*i-1,j+1)+rhh(2*i-2,j)+rhh(2*i,j))+1/4*(rhh(2*i-2,j+1)+rhh(2*i,j+1)));
            elseif i==m1_half %�ϱ߽�
                rH(i,j)=1/3*(rhh(2*i-1,2*j-1)+1/2*(rhh(2*i-2,2*j-1)+rhh(2*i-1,2*j-2)+rhh(2*i-1,2*j))+1/4*(rhh(2*i-2,2*j-2)+rhh(2*i-2,2*j)));
            elseif j==m2_half %�ұ߽�
                rH(i,j)=1/3*(rhh(2*i-1,2*j-1)+1/2*(rhh(2*i-2,2*j-1)+rhh(2*i-1,2*j-2)+rhh(2*i,2*j-1))+1/4*(rhh(2*i-2,2*j-2)+rhh(2*i,2*j-2)));
            else  %�ڵ�
                rH(i,j)=1/4*(rhh(2*i-1,2*j-1)+1/2*(rhh(2*i-2,2*j-1)+rhh(2*i-1,2*j-2)+rhh(2*i,2*j-1)+rhh(2*i-1,2*j))+1/4*(rhh(2*i-2,2*j-2)+rhh(2*i-2,2*j)+rhh(2*i,2*j-2)+rhh(2*i,2*j)));
            end
        end
    end
    %% �ڴ������Ͼ�ȷ����A1*eH=rH
    rHH=reshape(rH',N,1); %�Ѿ���ת��Ϊ����
    eH=A1\rHH;
    eHH=reshape(eH,m1_half,m2_half)';
    eh=zeros(m1-1,m2-1);
    %% ��У������ֵ��ϸ������
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
    %% �Բ�ֵ��Ľ��������й⻬����
    for i=1:3     % �ſ˱ȵ���3��
    u=-AD\(AL+AU)*u+AD\b;
    end
    %% solve_v
    u=reshape(u,m1-1,m2-1)';
    
    if( mod(lay,100) ==0 )  % ÿ��100��洢һ��
        fprintf('ʱ���Ϊ��%f \n',lay);
        eval(['save datau',num2str(lay)]);
    end    
    v = v + epsilon*tau*( beta * u -gamma * v -delta ); % ��ⳣ΢�ַ���
end
contour(u)