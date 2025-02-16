clear;clc;
df=xlsread('附件1.csv');
data=xlsread('附件2.csv');

%% 将所有坐标旋转
a=36.795*pi/180;        %方位角α
b=(90-78.169)*pi/180;   %仰角β
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

%% 旋转后求解
x1=df1(:,1);    %旋转后坐标
y1=df1(:,2);
z1=df1(:,3);
[col,row]=min(df1);
work=[];    %储存工作区域的主索节点
dir=[];     %储存相应的附件2的节点
bianhao=[]; %储存work中的节点对应的编号
k=1;
for i=1:2226
    if x1(i)^2+y1(i)^2<=150^2
        work(k,:)=[x1(i),y1(i),z1(i)];
        bianhao(k)=i;
        dir(k,:)=data1(i,:);
        k=k+1;
    end
end
bianhao=bianhao';
% [col,row]=max(dir); %边界点所在的位置
% for u=1:554
%     n1 = [dir(u,1)-dir(u,4);dir(u,2)-dir(u,5)]; %边界点促动器伸缩向量
%     n2 = [work(u,1);work(u,2)];   %边界点主索点坐标
%     gama(u) = acos(dot(n1,n2)/(norm(n1)*norm(n2))); %弧度制
%     gama(u) = gama(u)/pi*180; %换算成角度
% end
% gama=gama';
% gama=cos(gama);
t=0;
R=300.4;
m=1;
k=1;
ans_work=[];            %存储每个调节高度的主索节点可利用数量
ans_adjust=[];  %储存每个主索节点的调整量
ans_pnum=cell(1,15);    %储存调整后主索节点编号及坐标
for h=-0.6:0.1:0.6
% for h=-0.1:-0.1
    t=0;
    a = 1/((0.466*R-h+1.3313)*4);
    for i = 1:692
        x2=dir(i,1);
        y2=dir(i,2);
        z2=dir(i,3);
        x3=dir(i,4);
        y3=dir(i,5);
        z3=dir(i,6);
        syms x y z
        if x2==x3
            xx=0;
            [y,z]=solve((z2-z3)*(y-y2)/(y2-y3)+z2-z,a*y^2-R+h-z);
            if z(1)<0
                yy=y(1);
                zz=z(1);
            else
                yy=y(2);
                zz=z(2);
            end
            xx=double(xx);
            yy=double(yy);
            zz=double(zz);
            ans_pnum{k}(i,:)=[xx,yy,zz];
        else
            [x,y,z]=solve(a*(x^2+y^2)-R+h-z,(y2-y3)*(x-x2)/(x2-x3)+y2-y,(z2-z3)*(x-x2)/(x2-x3)+z2-z);
            if z(1)<0
                xx=x(1);
                yy=y(1);
                zz=z(1);
            else
                xx=x(2);
                yy=y(2);
                zz=z(2);
            end
            xx=double(xx);
            yy=double(yy);
            zz=double(zz);
            ans_pnum{k}(i,:)=[xx,yy,zz];
        end
        if xx^2+yy^2+zz^2<work(i,1)^2+work(i,2)^2+work(i,3)^2
            ans_adjust(i,k)=dis(xx,yy,zz,work(i,1),work(i,2),work(i,3));
            if ans_adjust(i,k)>0.6
                ans_adjust(i,k)=0.6;
            end
        else
            ans_adjust(i,k)=-dis(xx,yy,zz,work(i,1),work(i,2),work(i,3));
            if ans_adjust(i,k)<-0.6
                ans_adjust(i,k)=-0.6;
            end
        end
        if dis(xx,yy,zz,work(i,1),work(i,2),work(i,3))<=0.6
            t=t+1;
        end
    end
    ans_work(1,m)=t;
    ans_work(2,m)=h;
    m=m+1;
    k=k+1;
end
%绘制数据变化的可视化图
% load ans_work_p2.mat
% ans_work=ans_work_p2;
figure(2);
ax=ans_work(2,:);
ax(7)=0;
ay=ans_work(1,:)./692;
plot(ax,ay,'-');
text(-0.3,0.96,'(-0.2,0.9957)');
text(0,0.96,'(0,0.9957)');
xlabel('顶点的伸缩量h');
ylabel('主索节点的可行比例');
hold on;
plot(-0.2,0.9957,'r*');
plot(0,0.9957,'r*');
%% 基于以上寻优结果，继续在[-0.2,-0.01]处寻优
t=0;
R=300.4;
m=1;
k=1;
ans_work=[];            %存储每个调节高度的主索节点可利用数量
ans_adjust=[];  %储存每个主索节点的调整量
ans_pnum=cell(1,15);    %储存调整后主索节点编号及坐标
for h=-0.2:0.01:-0.01
    t=0;
    a = 1/((0.466*R-h+1.3313)*4);
    for i = 1:692
        x2=dir(i,1);
        y2=dir(i,2);
        z2=dir(i,3);
        x3=dir(i,4);
        y3=dir(i,5);
        z3=dir(i,6);
        syms x y z
        if x2==x3
            xx=0;
            [y,z]=solve((z2-z3)*(y-y2)/(y2-y3)+z2-z,a*y^2-R+h-z);
            if z(1)<0
                yy=y(1);
                zz=z(1);
            else
                yy=y(2);
                zz=z(2);
            end
            xx=double(xx);
            yy=double(yy);
            zz=double(zz);
            ans_pnum{k}(i,:)=[xx,yy,zz];
        else
            [x,y,z]=solve(a*(x^2+y^2)-R+h-z,(y2-y3)*(x-x2)/(x2-x3)+y2-y,(z2-z3)*(x-x2)/(x2-x3)+z2-z);
            if z(1)<0
                xx=x(1);
                yy=y(1);
                zz=z(1);
            else
                xx=x(2);
                yy=y(2);
                zz=z(2);
            end
            xx=double(xx);
            yy=double(yy);
            zz=double(zz);
            ans_pnum{k}(i,:)=[xx,yy,zz];
        end
        if xx^2+yy^2+zz^2<work(i,1)^2+work(i,2)^2+work(i,3)^2
            ans_adjust(i,k)=dis(xx,yy,zz,work(i,1),work(i,2),work(i,3));
            if ans_adjust(i,k)>0.6
                ans_adjust(i,k)=0.6;
            end
        else
            ans_adjust(i,k)=-dis(xx,yy,zz,work(i,1),work(i,2),work(i,3));
            if ans_adjust(i,k)<-0.6
                ans_adjust(i,k)=-0.6;
            end
        end
        if dis(xx,yy,zz,work(i,1),work(i,2),work(i,3))<=0.6
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
ay=ans_work(1,:)./692;
plot(ax,ay,'-');
xlabel('顶点的伸缩量h');
ylabel('主索节点的可行比例');

%% 将旋转后的坐标重置回原来位置
load ans_pnum.mat
a=36.795*pi/180;    %方位角α
b=(90-78.169)*pi/180;    %仰角β
Ry_ni=[cos(b) 0 -sin(b)
    0 1 0
    sin(b) 0 cos(b)];
Rz_ni=[cos(a) sin(a) 0
    -sin(a) cos(a) 0
    0 0 1];
for i = 1:692
    res(i,:)=ans_pnum{6}(i,:)*Ry_ni*Rz_ni;
end
plot3(res(:,1),res(:,2),res(:,3),'*');
xlabel('x');
ylabel('y');
zlabel('z');
%中心主索节点
[col,row]=min(ans_pnum{6});
res(132,:)

%% 附件1所有主索节点的编号
A = readmatrix('附件1.csv', 'OutputType', 'string');
A(:,2:4)=[];        %储存所有节点编号
B=strings(692,1);   %储存相应的节点编号
for i=1:692
    B(i)=A(bianhao(i));
end

%% 计算旋转后的调节量
% for i=1:692
%     jizhun(i,:)=work(i,:)*Ry_ni*Rz_ni;
% end
% d=[];   %储存调节量
% for i=1:692
%     if res(i,1)^2+res(i,2)^2+res(i,3)^2<jizhun(i,1)^2+jizhun(i,2)^2+jizhun(i,3)
%         d(i)=dis(res(i,1),res(i,2),res(i,3),jizhun(i,1),jizhun(i,2),jizhun(i,3));
%         if d(i)>0.6
%             d(i)=0.6;
%         end
%     else
%         d(i)=-dis(res(i,1),res(i,2),res(i,3),jizhun(i,1),jizhun(i,2),jizhun(i,3));
%         if d(i)<-0.6
%             d(i)=-0.6;
%         end
%     end
% end
% d=d';

%% 计算所有有用主索节点对应的反射板数量
load ans_adjust.mat
res_adjust=ans_adjust(:,6);
for i=692:-1:1  %从后往前遍历，防止数据顶替
    if res_adjust(i)>0.6 || res_adjust(i)<-0.6
        B(i)=[];
    end
end

dt = readmatrix('附件3.csv', 'OutputType', 'string');
cishu=0;    
bankuai=0;  %储存板块数
for i=2:4301
    for j=1:3
        for k=1:692
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




