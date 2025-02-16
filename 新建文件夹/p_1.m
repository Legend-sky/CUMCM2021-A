clear;clc;
df=xlsread('附件1.csv');
data=xlsread('附件2.csv');
x=df(:,1);
y=df(:,2);
z=df(:,3);
% createFit1(x, y, z); %拟合图像
figure(1);
plot3(x,y,z,'*','Color',[0 1 1]);
hold on;
plot3(data(:,1),data(:,2),data(:,3),'*');
hold on;
plot3(data(:,4),data(:,5),data(:,6),'*');
hold on;
grid on;
%% 最优抛物面
h=0;    %首先绘制出未调整时的最优抛物面
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
legend('基准球面','促动器下端点位置','基准态上端点');
%% 寻找工作区域
work=[];    %储存工作区的节点坐标
dir=[];     %储存相应附件2的信息
bianhao=[]; %储存work中的节点对应的编号
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
%% 促动器向量与径向夹角
[col,row]=max(dir); %边界点所在的位置
gama(1)=0;
for u=2:706
    n1 = [dir(u,1)-dir(u,4);dir(u,2)-dir(u,5)]; %边界点促动器伸缩向量
    n2 = [work(u,1);work(u,2)];   %边界点主索点坐标
    gama(u) = acos(dot(n1,n2)/(norm(n1)*norm(n2))); %弧度制
    gama(u) = gama(u)/pi*180; %换算成角度
end
gama=gama';
gama=cos(gama);
%% 计算最优抛物面
% k = work(row(3),3)/work(row(3),1);
% outx=work(row(3),1);    %边界节点的x坐标
% outy=work(row(3),2);    %边界节点的y坐标
% syms px;
% f=k*px-4*(0.466*R-h)*px^2+300.4-h;
% px=vpa(solve(f==0),4)
k=1;
t=0;
m=1;
ans_work=[];    %存储每个调节高度的主索节点可利用数量
ans_adjust=[];  %储存每个主索节点的调整量
ans_pnum=cell(1,15); %储存调整后的住主索节点坐标
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
            wy=solve(a.*wy.^2-z2./y2.*wy-R+h);  %求解直线与抛物面相交的交点
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
%绘制数据变化的可视化图
figure(2);
ax=ans_work(2,:);
ax(7)=0;
ay=ans_work(1,:)./706;
plot(ax,ay,'-');
xlabel('调整的h高度');
ylabel('主索节点的可行比例');

%% 在以上基础上迭代寻优，范围是[-0.3，-0.01]
k=1;
t=0;
m=1;
ans_work=[];    %存储每个调节高度的主索节点可利用数量
ans_adjust=[];  %储存每个主索节点的调整量

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
%绘制数据变化的可视化图
figure(3);
ax=ans_work(2,:);
ay=ans_work(1,:)./706;
plot(ax,ay,'-');
xlabel('调整的h高度');
ylabel('主索节点的可行比例');

%% 计算对应的反射面板数量
A = readmatrix('附件1.csv', 'OutputType', 'string');
A(:,2:4)=[];        %储存所有节点编号
B=strings(706,1);   %储存相应的节点编号
for i=1:706
    B(i)=A(bianhao(i));
end

res_adjust=ans_adjust(:,6);
for i=706:-1:1  %从后往前遍历，防止数据顶替
    if res_adjust(i)>0.6 || res_adjust(i)<-0.6
        B(i)=[];
    end
end
dt = readmatrix('附件3.csv', 'OutputType', 'string');
cishu=0;    
bankuai=0;  %储存板块数
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

%% 最优抛物面与工作点的相交图
h=-0.15;
R=300.4;
% z1=(0.466*R-h)*4.*(x1.^2+y1.^2)-300.4+h;
lixiang=[];     %理想抛物面
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
legend('基准球面','工作抛物面');
grid on;

