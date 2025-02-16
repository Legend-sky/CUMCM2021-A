clear;clc;
df=xlsread('附件1.csv');
data=xlsread('附件2.csv');
% 将所有坐标旋转
a=36.795*pi/180;        %方位角α
b=(90-78.169)*pi/180;   %仰角β
R=300.4;
Ry=[cos(b) 0 sin(b)
    0 1 0
    -sin(b) 0 cos(b)];
Rz=[cos(a) -sin(a) 0
    sin(a) cos(a) 0
    0 0 1];
Ry_ni=[cos(b) 0 -sin(b)
    0 1 0
    sin(b) 0 cos(b)];
Rz_ni=[cos(a) sin(a) 0
    -sin(a) cos(a) 0
    0 0 1];
for i = 1:2226
    df1(i,:)=df(i,:)*Rz*Ry;
    data1(i,1:3)=data(i,1:3)*Rz*Ry;
    data1(i,4:6)=data(i,4:6)*Rz*Ry;
end
A = readmatrix('附件1.csv', 'OutputType', 'string');
A(:,2:4)=[];        %储存所有节点编号
dt = readmatrix('附件3.csv', 'OutputType', 'string');
w=1;
for D=250:10:350
    % 旋转后求解
    x1=df1(:,1);    %旋转后坐标
    y1=df1(:,2);
    z1=df1(:,3);
    [col,row]=min(df1);
    work=[];    %储存工作区域的主索节点
    dir=[];     %储存相应的附件2的节点
    bianhao=[]; %储存work中的节点对应的编号
    k=1;
    for i=1:2226
        if x1(i)^2+y1(i)^2<=(D/2)^2
            work(k,:)=[x1(i),y1(i),z1(i)];
            bianhao(k)=i;
            dir(k,:)=data1(i,:);
            k=k+1;
        end
    end
    bianhao=bianhao';
    n=size(bianhao,1);
    t=0;
    m=1;
    k=1;
    ans_work=[];            %存储每个调节高度的主索节点可利用数量
    ans_adjust=[];          %储存每个主索节点的调整量
    ans_pnum=cell(1,15);    %储存调整后主索节点编号及坐标
    for h=-0.6:0.1:0.6
        t=0;
        a = 1/((0.466*R-h+1.3313)*4);
        for i = 1:n
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
            else
                ans_adjust(i,k)=-dis(xx,yy,zz,work(i,1),work(i,2),work(i,3));
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
    [val,index]=max(ans_work(1,:)); %寻找其中最优的点，在index列
    
    S_zong=0;

    % 读取调整后节点的坐标
    for i = 1:n
        res(i,:)=ans_pnum{index}(i,:)*Ry_ni*Rz_ni;  %储存调节之后的节点坐标
    end
    
    B=strings(n,1);   %储存相应的节点编号
    for i=1:n
        B(i)=A(bianhao(i));
    end
    % 计算总面积
    
    cishu=0;    
    bankuai=0;  %储存板块数
    arf=36.795*pi/180;
    beta=78.169*pi/180;
    n2=[cos(beta)*cos(arf),cos(beta)*sin(arf),sin(beta)];
    u=1;
    for i=2:4301
        for j=1:3
            for k=1:n
                if dt(i,j)==B(k)
                    cishu=cishu+1;
                    break;
                end
            end
            if j==1
                a=res(k,:);
            elseif j==2
                b=res(k,:);
            else
                c=res(k,:);
            end
        end
        if cishu==3
            a_temp=(b(2)-a(2))*(c(3)-a(3))-(b(3)-a(3))*(c(2)-a(2));
            b_temp=(b(3)-a(3))*(c(1)-a(1))-(b(1)-a(1))*(c(3)-a(3));
            c_temp=(b(1)-a(1))*(c(2)-a(2))-(b(2)-a(2))*(c(1)-a(1));
            n1=[a_temp,b_temp,c_temp];
            gama(u) = dot(n1,n2)/(norm(n1)*norm(n2)); %弧度制
            S_zong=S_zong+Area(a,b,c)*abs(gama(u));
            u=u+1;
        end
        cishu=0;
    end
    
    arf=36.795*pi/180;
    beta=78.169*pi/180;
    n2=[cos(beta)*cos(arf),cos(beta)*sin(arf),sin(beta)];
    res_adjust=ans_adjust(:,index);
    B3=B;
    res3=res;
    kk=0;
    for i=n:-1:1  %从后往前遍历，防止数据顶替
        if res_adjust(i)>0.6 || res_adjust(i)<-0.6
            B3(i)=[];
            res3(i,:)=[];
            kk=kk+1;
        end
    end
    u=1;
    S_work=0;
    for i=2:4301
        for j=1:3
            for k=1:n-kk
                if dt(i,j)==B3(k)
                    cishu=cishu+1;
                    break;
                end
            end
            if j==1
                a=res3(k,:);
            elseif j==2
                b=res3(k,:);
            else
                c=res3(k,:);
            end
        end
        if cishu==3
            a_temp=(b(2)-a(2))*(c(3)-a(3))-(b(3)-a(3))*(c(2)-a(2));
            b_temp=(b(3)-a(3))*(c(1)-a(1))-(b(1)-a(1))*(c(3)-a(3));
            c_temp=(b(1)-a(1))*(c(2)-a(2))-(b(2)-a(2))*(c(1)-a(1));
            n1=[a_temp,b_temp,c_temp];
            gama3(u) = dot(n1,n2)/(norm(n1)*norm(n2)); %弧度制
            S_work=S_work+Area(a,b,c)*abs(gama3(u));
            u=u+1;
        end
        cishu=0;
    end
    use=S_work/S_zong;
    jieguo(w,:)=[D,S_work,S_zong,use];
    w=w+1;
end
plot(jieguo(:,1),jieguo(:,4),'r-');


