function NSFDM_MG()
m1 = 256;
m2 = 256;
m1_half = m1/2;
m2_half = m2/2;
N = (m1_half-1)*(m2_half-1);
M = (m1 - 1)*(m2 - 1);
a = 0.1;
alpha1 = 2;
alpha2 = 2;
Dx = 2.5;
Dy = 2.5;
tau = 0.1;
Kx = 0.0001;
Ky = 0.0001;
epsilon = 0.01;
beta = 0.5;
gamma = 1;
delta = 0;
t = 1000;
b = zeros(m1-1,m2-1);
u = zeros(m1-1,m2-1);
v = zeros(m1-1,m2-1);
g1 = zeros(1,m1+1);
g2 = zeros(1,m1+1);
hx = Dx/m1; %ϸ���񲽳�
hy = Dy/m2; 
Hx = Dx/m1_half; %�����񲽳�
Hy = Dy/m2_half;

Phi1t=(1-exp(-a*tau))/a;

Ca1 = 1.0/(2 * cos(pi * alpha1 / 2));
Ca2 = 1.0/(2 * cos(pi * alpha2 / 2));

r1 = (-1) * (Phi1t* Kx * Ca1)/(hx^alpha1);
r2 = (-1) * (Phi1t * Ky * Ca2)/(hy^alpha2);
rr1= (-1) * (Phi1t * Kx * Ca1)/(Hx^alpha1);
rr2= (-1) * (Phi1t* Ky * Ca2)/(Hy^alpha2);

layend = t / tau;
xx=linspace(0,Dx,m1-1);
yy=linspace(0,Dy,m2-1);
xx=repmat(xx,1,m1-1)';
yy=repmat(yy,m2-1,1);
yy=reshape(yy,M,1);
%% initial u,v 
for i = 1:m1-1
    for j = 1:m2-1
        if( ( i < m1_half ) && ( j <= m2_half ) )
            u( i,j ) = 1.0;
        end
        if( ( i >= m1_half ) && ( i < m1 ) )
            v( i,j  ) = 0.1;
        end
    end
end
%% initial g(L)
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
 for i = 1:m1_half-1
    for j = 1:m2_half-1
        A1((i-1)*(m2_half-1)+j,(i-1)*(m2_half-1)+j) = 1.0; %ϵ������Խ����ϵ�ֵ
        if i == m1_half-1    % ��������߽����⴦��
                A1((i-1)*(m2_half-1)+j,(i-(1:i)+1-1)*(m2_half-1)+j) = A1((i-1)*(m2_half-1)+j,(i-(1:i)+1-1)*(m2_half-1)+j)+(-1*rr1 * g1((1:i)+1));
        else
                A1((i-1)*(m2_half-1)+j,(i-(0:i)+1-1)*(m2_half-1)+j) = A1((i-1)*(m2_half-1)+j,(i-(0:i)+1-1)*(m2_half-1)+j)+(-1*rr1 * g1((0:i)+1));
        end
        if i == 1          % �������ұ߽����⴦��
                A1((i-1)*(m2_half-1)+j,(i+(1:m1_half-1-i)-1-1)*(m2_half-1)+j) = A1((i-1)*(m2_half-1)+j,(i+(1:m1_half-1-i)-1-1)*(m2_half-1)+j)+(-1*rr1 * g1((1:m1_half-1-i)+1));
        else
                A1((i-1)*(m2_half-1)+j,(i+(0:m1_half-1-i)-1-1)*(m2_half-1)+j) = A1((i-1)*(m2_half-1)+j,(i+(0:(m2_half-1)-i)-1-1)*(m2_half-1)+j)+(-1*rr1 * g1((0:m1_half-1-i)+1));
        end
        if j == m2_half-1    % �������±߽����⴦��
                A1((i-1)*(m2_half-1)+j,(i-1)*(m2_half-1)+j-(1:j)+1) = A1((i-1)*(m2_half-1)+j,(i-1)*(m2_half-1)+j-(1:j)+1)+(-1*rr2 * g2((1:j)+1));
        else
                A1((i-1)*(m2_half-1)+j,(i-1)*(m2_half-1)+j-(0:j)+1) = A1((i-1)*(m2_half-1)+j,(i-1)*(m2_half-1)+j-(0:j)+1)+(-1*rr2 * g2((0:j)+1));
        end
        if j == 1          % �������ϱ߽����⴦��
                A1((i-1)*(m2_half-1)+j,(i-1)*(m2_half-1)+j+(1:m1_half-1-j)-1) = A1((i-1)*(m2_half-1)+j,(i-1)*(m2_half-1)+j+(1:m1_half-1-j)-1)+(-1*rr2 * g2((1:m1_half-1-j)+1));
        else
                A1((i-1)*(m2_half-1)+j,(i-1)*(m2_half-1)+j+(0:m1_half-1-j)-1) = A1((i-1)*(m2_half-1)+j,(i-1)*(m2_half-1)+j+(0:m1_half-1-j)-1)+(-1*rr2 * g2((0:m1_half-1-j)+1));
        end
    end
    fprintf('��ɾ���A1��������%f \n',i);     % ÿ��ɾ���A��һ�����һ��
 end
fprintf('Matrix A1 is over\n');      % ����A�������
save A1 A1;


for lay = 1:layend
     % comute the right-hand side b
     b=u+2*Phi1t*u.*u.*u-Phi1t*(1+a)*u.*u-Phi1t*v;
       
     % update the fine-grid matrix A
     Adiag=3*Phi1t*u.*u-2*Phi1t*(1+a)*u+a*Phi1t;
     A=A+sparse(1:M,1:M,Adiag'); 
 
     % update the coarse-grid matrix A1
     A1diag=(1+a)*Phi1t*u+Phi1t*u.*u+2*a*Phi1t;
     A1=A1+sparse(1:N,1:N,A1diag(2:2:end,2:2:end));
    
    b=reshape(b',M,1);
    u=reshape(u',M,1); %��u�����ų�������
    AL=tril(A,-1);%ϵ������A����������
    AU=triu(A,1);%ϵ������A����������
    AD=diag(diag(A));
    for i=1:3
    u=-AD\(AL+AU)*u+AD\b;
    end
    
    %% ��������rh
    rh=b-A*u;
    %% ��������������������ΪrH
    rhh=reshape(rh,m1-1,m2-1)';
    rH=zeros(m1_half-1,m2_half-1);
    for i=1:m1_half-1
        for j=1:m2_half-1
            if i==1&&j==1 %���½ǵ�
                rH(i,j)=4/9*(rhh(i,j)+1/2*(rhh(i,j+1)+rhh(i+1,j))+1/4*rhh(i+1,j+1));
            elseif i==1&&j==m2_half-1  %���½ǵ�
                rH(i,j)=4/9*(rhh(i,2*j-1)+1/2*(rhh(i,2*j-2)+rhh(i+1,2*j-1))+1/4*rhh(i+1,2*j-2));
            elseif i==m1_half-1&&j==1  %���Ͻǵ�
                rH(i,j)=4/9*(rhh(2*i-1,j)+1/2*(rhh(2*i-2,j)+rhh(2*i-1,j+1))+1/4*rhh(2*i-2,j+1));
            elseif i==m1_half-1&&j==m2_half-1 %���Ͻǵ�
                 rH(i,j)=4/9*(rhh(2*i-1,2*j-1)+1/2*(rhh(2*i-1,2*j-2)+rhh(2*i-2,2*j-1))+1/4*rhh(2*i-2,2*j-2));
            elseif i==1 %�±߽�
                rH(i,j)=1/3*(rhh(i,2*j-1)+1/2*(rhh(i,2*j-2)+rhh(i,2*j)+rhh(i+1,2*j-1))+1/4*(rhh(i+1,2*j-2)+rhh(i+1,2*j)));
            elseif j==1 %��߽�
                rH(i,j)=1/3*(rhh(2*i-1,j)+1/2*(rhh(2*i-1,j+1)+rhh(2*i-2,j)+rhh(2*i,j))+1/4*(rhh(2*i-2,j+1)+rhh(2*i,j+1)));
            elseif i==m1_half-1 %�ϱ߽�
                rH(i,j)=1/3*(rhh(2*i-1,2*j-1)+1/2*(rhh(2*i-2,2*j-1)+rhh(2*i-1,2*j-2)+rhh(2*i-1,2*j))+1/4*(rhh(2*i-2,2*j-2)+rhh(2*i-2,2*j)));
            elseif j==m2_half-1 %�ұ߽�
                rH(i,j)=1/3*(rhh(2*i-1,2*j-1)+1/2*(rhh(2*i-2,2*j-1)+rhh(2*i-1,2*j-2)+rhh(2*i,2*j-1))+1/4*(rhh(2*i-2,2*j-2)+rhh(2*i,2*j-2)));
            else  %�ڵ�
                rH(i,j)=1/4*(rhh(2*i-1,2*j-1)+1/2*(rhh(2*i-2,2*j-1)+rhh(2*i-1,2*j-2)+rhh(2*i,2*j-1)+rhh(2*i-1,2*j))+1/4*(rhh(2*i-2,2*j-2)+rhh(2*i-2,2*j)+rhh(2*i,2*j-2)+rhh(2*i,2*j)));
            end
        end
    end
    %% �ڴ������Ͼ�ȷ����A1*eH=rH
    rHH=reshape(rH',N,1); %�Ѿ���ת��Ϊ����
    eH=A1\rHH;
    eHH=reshape(eH,m1_half-1,m2_half-1)';
    eh=zeros(m1-1,m2-1);
    %��У������ֵ��ϸ������
    for i=1:m1_half-1
        for j=1:m2_half-1
            if i<m1_half-1&&j<m2_half-1
                eh(2*i-1,2*j-1)=eHH(i,j);
                eh(2*i-1,2*j)=1/2*(eHH(i,j)+eHH(i,j+1));
                eh(2*i,2*j-1)=1/2*(eHH(i,j)+eHH(i+1,j));  
                eh(2*i,2*j)=1/4*(eHH(i,j)+eHH(i,j+1)+eHH(i+1,j)+eHH(i+1,j+1));
            elseif i==m1_half-1&&j<m2_half-1
                    eh(2*i-1,2*j-1)=eHH(i,j);
                    eh(2*i-1,2*j)=1/2*(eHH(i,j)+eHH(i,j+1));
            elseif j==m2_half-1&&i<m1_half-1
                    eh(2*i-1,2*j-1)=eHH(i,j);
                    eh(2*i,2*j-1)=1/2*(eHH(i,j)+eHH(i+1,j)); 
            end
        end
    end
    ehh=reshape(eh',M,1);
    u=u+ehh;
    %�Բ�ֵ��Ľ��������й⻬�������ſ˱ȵ���3��
    for i=1:3
    u0=-AD\(AL+AU)*u+AD\b;
    end
    %%solve_v
    u=reshape(u0,m1-1,m2-1)';
    v = v + epsilon*tau*( beta * u -gamma * v -delta );
    if( mod(lay,100) ==0 )
%         uu=[xx,yy,u0];
%         [row,col]=size(uu);
%         str=['data_',num2str(lay),'.plt'];
%         fid=fopen(str,'w+');
%         fprintf(fid,'TITLE="contour of the solution"\n');
%         fprintf(fid,'VARIABLES="x","y","u"\n');
%         fprintf(fid,'ZONE T="BOX",I= 255,J= 255,F=POINT\n');
%         for i=1:row
%             for j=1:col
%                 if j==col
%                     fprintf(fid,'%g\n',uu(i,j));
%                 else
%                     fprintf(fid,'%g\t',uu(i,j));
%                 end
%             end
%         end
%         fclose(fid);
        fprintf('ʱ���Ϊ��%f \n',lay);
    end
end
contour(u)