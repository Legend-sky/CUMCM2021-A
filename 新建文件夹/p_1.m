clear;clc;
df=xlsread('����1.csv');
data=xlsread('����2.csv');
x=df(:,1);
y=df(:,2);
z=df(:,3);
% createFit1(x, y, z); %���ͼ��
figure(1);
plot3(x,y,z,'*','Color',[0 1 1]);
hold on;
plot3(data(:,1),data(:,2),data(:,3),'*');
hold on;
plot3(data(:,4),data(:,5),data(:,6),'*');
hold on;
grid on;
%% ����������
h=0;    %���Ȼ��Ƴ�δ����ʱ������������
R=300.4;
x1=[-250:1:250]';
y1=[-250:1:250]';
% z1=(0.466*R-h)*4.*(x1.^2+y1.^2)-300.4+h;
z1=[];
for i=1:501
    for j=1:501
        z1(i,j)=(1/((0.466*R-h)*4)).*(x1(i).^2+y1(j).^2)-300.4+h;
    end
end
mesh(x1,y1,z1);
xlabel('x');
ylabel('y');
zlabel('z');
legend('��׼����','�ٶ����¶˵�λ��','��׼̬�϶˵�');
%% Ѱ�ҹ�������
work=[];    %���湤�����Ľڵ�����
dir=[];     %������Ӧ����2����Ϣ
bianhao=[]; %����work�еĽڵ��Ӧ�ı��
k=1;
for i=1:2226
    if x(i)^2+y(i)^2 <= 150^2
        work(k,:)=[x(i),y(i),z(i)];
        bianhao(k)=i;
        dir(k,:)=data(i,:);
        k=k+1;
    end
end
bianhao=bianhao';
%% �ٶ��������뾶��н�
[col,row]=max(dir); %�߽�����ڵ�λ��
gama(1)=0;
for u=2:706
    n1 = [dir(u,1)-dir(u,4);dir(u,2)-dir(u,5)]; %�߽��ٶ�����������
    n2 = [work(u,1);work(u,2)];   %�߽������������
    gama(u) = acos(dot(n1,n2)/(norm(n1)*norm(n2))); %������
    gama(u) = gama(u)/pi*180; %����ɽǶ�
end
gama=gama';
gama=cos(gama);
%% ��������������
% k = work(row(3),3)/work(row(3),1);
% outx=work(row(3),1);    %�߽�ڵ��x����
% outy=work(row(3),2);    %�߽�ڵ��y����
% syms px;
% f=k*px-4*(0.466*R-h)*px^2+300.4-h;
% px=vpa(solve(f==0),4)
k=1;
t=0;
m=1;
ans_work=[];    %�洢ÿ�����ڸ߶ȵ������ڵ����������
ans_adjust=[];  %����ÿ�������ڵ�ĵ�����
ans_pnum=cell(1,15); %����������ס�����ڵ�����
for h=-0.6:0.1:0.6
% for h=-0.15:-0.15
    t=0;
    a = 1/((0.466*R-h+1.3313)*4);
    ans_adjust(1,k)=h;
    ans_pnum{k}(1,:)=[0,0,-300.4+h];
    for i = 2:706
        x2=work(i,1);
        y2=work(i,2);
        z2=work(i,3);
        syms wx wy wz
        if x2==0
            wy=solve(a.*wy.^2-z2./y2.*wy-R+h);  %���ֱ�����������ཻ�Ľ���
            wz=z2./y2.*wy;
            if wz(1)<0
                yy=wy(1);
                zz=wz(1);
            else
                yy=wy(2);
                zz=wz(2);
            end
            xx=double(xx);
            yy=double(yy);
            zz=double(zz);
            ans_pnum{k}(i,:)=[xx,yy,zz];
        else
            wx=solve(a.*(1+y2.^2./x2.^2).*wx.^2-z2./x2.*wx-R+h);        
            wy=y2./x2.*wx;
            wz=z2./x2.*wx;
            if wz(1)<0
                xx=wx(1);
                yy=wy(1);
                zz=wz(1);
            else
                xx=wx(2);
                yy=wy(2);
                zz=wz(2);
            end
            xx=double(xx);
            yy=double(yy);
            zz=double(zz);
            ans_pnum{k}(i,:)=[xx,yy,zz];
        end
        if xx^2+yy^2+zz^2<x2^2+y2^2+z2^2
            ans_adjust(i,k)=dis(xx,yy,zz,x2,y2,z2);
            if ans_adjust(i,k)>0.6
                ans_adjust(i,k)=0.6;
            end
        else
            ans_adjust(i,k)=-dis(xx,yy,zz,x2,y2,z2);
            if ans_adjust(i,k)<-0.6
                ans_adjust(i,k)=-0.6;
            end
        end
        if dis(xx,yy,zz,x2,y2,z2)<=0.6
            t=t+1;
        end
    end
    ans_work(1,m)=t;
    ans_work(2,m)=h;
    m=m+1;
    k=k+1;
end
%�������ݱ仯�Ŀ��ӻ�ͼ
figure(2);
ax=ans_work(2,:);
ax(7)=0;
ay=ans_work(1,:)./706;
plot(ax,ay,'-');
xlabel('������h�߶�');
ylabel('�����ڵ�Ŀ��б���');

%% �����ϻ����ϵ���Ѱ�ţ���Χ��[-0.3��-0.01]
k=1;
t=0;
m=1;
ans_work=[];    %�洢ÿ�����ڸ߶ȵ������ڵ����������
ans_adjust=[];  %����ÿ�������ڵ�ĵ�����

for h=-0.3:0.02:0.01
    t=0;
    a = 1/((0.466*R-h+1.3313)*4);
    ans_pnum{k}(1,:)=[0,0,-300.4+h];
    for i = 2:706
        x2=work(i,1);
        y2=work(i,2);
        z2=work(i,3);
        syms wx wy wz
        if x2==0
            wy=solve(a.*wy.^2-z2./y2.*wy-R+h);
            wz=z2./y2.*wy;
            if wz(1)<0
                yy=wy(1);
                zz=wz(1);
            else
                yy=wy(2);
                zz=wz(2);
            end
            xx=double(xx);
            yy=double(yy);
            zz=double(zz);
            ans_pnum{k}(i,:)=[xx,yy,zz];
        else
            wx=solve(a.*(1+y2.^2./x2.^2).*wx.^2-z2./x2.*wx-R+h);        
            wy=y2./x2.*wx;
            wz=z2./x2.*wx;
            if wz(1)<0
                xx=wx(1);
                yy=wy(1);
                zz=wz(1);
            else
                xx=wx(2);
                yy=wy(2);
                zz=wz(2);
            end
            xx=double(xx);
            yy=double(yy);
            zz=double(zz);
            ans_pnum{k}(i,:)=[xx,yy,zz];
        end
        if xx^2+yy^2+zz^2<x2^2+y2^2+z2^2
            ans_adjust(i,k)=dis(xx,yy,zz,x2,y2,z2);
            if ans_adjust(i,k)>0.6
                ans_adjust(i,k)=0.6;
            end
        else
            ans_adjust(i,k)=-dis(xx,yy,zz,x2,y2,z2);
            if ans_adjust(i,k)<-0.6
                ans_adjust(i,k)=-0.6;
            end
        end
        if dis(xx,yy,zz,x2,y2,z2)<=0.6
            t=t+1;
        end
    end
    ans_work(1,m)=t;
    ans_work(2,m)=h;
    m=m+1;
    k=k+1;
end
%�������ݱ仯�Ŀ��ӻ�ͼ
figure(3);
ax=ans_work(2,:);
ay=ans_work(1,:)./706;
plot(ax,ay,'-');
xlabel('������h�߶�');
ylabel('�����ڵ�Ŀ��б���');

%% �����Ӧ�ķ����������
A = readmatrix('����1.csv', 'OutputType', 'string');
A(:,2:4)=[];        %�������нڵ���
B=strings(706,1);   %������Ӧ�Ľڵ���
for i=1:706
    B(i)=A(bianhao(i));
end

res_adjust=ans_adjust(:,6);
for i=706:-1:1  %�Ӻ���ǰ��������ֹ���ݶ���
    if res_adjust(i)>0.6 || res_adjust(i)<-0.6
        B(i)=[];
    end
end
dt = readmatrix('����3.csv', 'OutputType', 'string');
cishu=0;    
bankuai=0;  %��������
for i=2:4301
    for j=1:3
        for k=1:683
            if dt(i,j)==B(k)
                cishu=cishu+1;
                break;
            end
        end
    end
    if cishu==3
        bankuai=bankuai+1;
    end
    cishu=0;
end

%% �����������빤������ཻͼ
h=-0.15;
R=300.4;
% z1=(0.466*R-h)*4.*(x1.^2+y1.^2)-300.4+h;
lixiang=[];     %����������
for i=1:706
    lixiang(i,1)=work(i,1);
    lixiang(i,2)=work(i,2);
    lixiang(i,3)=(1/((0.466*R-h)*4)).*(lixiang(i,1).^2+lixiang(i,2).^2)-300.4+h;
end

plot3(work(:,1),work(:,2),work(:,3),'*','Color',[0 1 1]);
hold on;
xlabel('x');
ylabel('y');
zlabel('z');
plot3(lixiang(:,1),lixiang(:,2),lixiang(:,3),'*','Color',[1	0.5 0]);
legend('��׼����','����������');
grid on;

