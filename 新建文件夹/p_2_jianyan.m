clear;clc;
df=xlsread('附件1.csv');
data=xlsread('附件2.csv');

a=0*pi/180;        %方位角α
b=(90-90)*pi/180;   %仰角β
Ry=[cos(b) 0 sin(b)
    0 1 0
    -sin(b) 0 cos(b)];
Rz=[cos(a) -sin(a) 0
    sin(a) cos(a) 0
    0 0 1];
for i = 1:2226
    df1(i,:)=df(i,:)*Rz*Ry;
    data1(i,1:3)=data(i,1:3)*Rz*Ry;
    data1(i,4:6)=data(i,4:6)*Rz*Ry;
end
figure(1);
plot3(df1(:,1),df1(:,2),df1(:,3),'*');
hold on;
plot3(data1(:,1),data1(:,2),data1(:,3));
hold on;
plot3(data1(:,4),data1(:,5),data1(:,6));
grid on;
xlabel('x');
ylabel('y');
zlabel('z');

work=[];    %储存工作区的节点坐标
dir=[];     %储存相应附件2的信息
k=1;
for i=1:2226
    if df1(i,1)^2+df1(i,2)^2 <= 150^2
        work(k,:)=[df1(i,1),df1(i,2),df1(i,3)];
        dir(k,:)=data1(i,:);
        k=k+1;
    end
end
R=300.4;
t=0;
m=1;
ans_work=[];    %存储每个调节高度的主索节点可利用数量
for h=-0.6:0.1:0.6
    t=0;
    a = 1/((0.466*R-h+1.3313)*4);
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
        end
        if dis(xx,yy,zz,x2,y2,z2)<=0.6
            t=t+1;
        end
    end
    ans_work(1,m)=t;
    ans_work(2,m)=h;
    m=m+1;
end
%绘制数据变化的可视化图
figure(2);
ax=ans_work(2,:);
ax(7)=0;
ay=ans_work(1,:)./706;
plot(ax,ay,'-');
xlabel('顶点的伸缩量h');
ylabel('主索节点的可行比例');

