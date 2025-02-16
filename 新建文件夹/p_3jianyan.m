%% 利用问题一的数据求解接受比
clear;clc;
load ans_adjust_p1.mat
load ans_pnum_p1.mat
load bianhao_p2.mat

S_zong=0;

for i = 1:706
    res(i,:)=ans_pnum{1}(i,:);  %储存调节之后的节点坐标
end

A = readmatrix('附件1.csv', 'OutputType', 'string');
A(:,2:4)=[];        %储存所有节点编号
B=strings(706,1);   %储存相应的节点编号
for i=1:706
    B(i)=A(bianhao(i));
end
% 计算总面积
dt = readmatrix('附件3.csv', 'OutputType', 'string');
cishu=0;    
bankuai=0;  %储存板块数
n2=[0,0,1];
u=1;
for i=2:4301
    for j=1:3
        for k=1:706
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
        vec(u,:)=n1;
        gama(u) = dot(n1,n2)/(norm(n1)*norm(n2)); %弧度制
        S_zong=S_zong+Area(a,b,c)*abs(gama(u));
        area1(u)=Area(a,b,c);
        u=u+1;
    end
    cishu=0;
end
area1=area1';
gama=gama';

%% 基准球面面积
S_jizhun=0;     %基准球面面积
R=300.4;
t=1;
min_ans=0;
f=1;
for i=-298:0.1:-2
    xita1=atan((-160.4136+sqrt(R^2-i^2))/(-0.5-i))*180/pi;
    n1=[-0.5-i,-160.4136+sqrt(R^2-i^2)];
    k=-sqrt(R^2-i^2)/i;
    n2=[1/sqrt(1+k^2),k/sqrt(1+k^2)];
    xita2=acos(dot(n1,n2)/(norm(n1)*norm(n2)))*180/pi;
    mk(t,1)=i;
    mk(t,2)=xita1+2*xita2;
    if f && mk(t,2)>=90
        f=0;
        min_ans=mk(t,1);
    end
    t=t+1;
end
min_ans=abs(min_ans-4);     %考虑到一块板的长度大约是4米
bankuai2=0;     %有效反射面板数量
res2=[];        %储存基准面上的有效主索节点
B2=strings(351,1);
k=1;
for i = 1:706
    if res(i,1)^2+res(i,2)^2<=min_ans^2
        res2(k,:)=res(i,:);
        B2(k)=B(i);
        k=k+1;
    end
end
u=1;
n2=[0,0,1];
for i=2:4301
    for j=1:3
        for k=1:351
            if dt(i,j)==B2(k)
                cishu=cishu+1;
                break;
            end
        end
        if j==1
            a=res2(k,:);
        elseif j==2
            b=res2(k,:);
        else
            c=res2(k,:);
        end
    end
    if cishu==3
        a_temp=(b(2)-a(2))*(c(3)-a(3))-(b(3)-a(3))*(c(2)-a(2));
        b_temp=(b(3)-a(3))*(c(1)-a(1))-(b(1)-a(1))*(c(3)-a(3));
        c_temp=(b(1)-a(1))*(c(2)-a(2))-(b(2)-a(2))*(c(1)-a(1));
        n1=[a_temp,b_temp,c_temp];
        vec2(u,:)=n1;
        gama2(u) = dot(n1,n2)/(norm(n1)*norm(n2)); %弧度制
        S_jizhun=S_jizhun+Area(a,b,c)*abs(gama2(u));
        area2(u)=Area(a,b,c);
        u=u+1;
    end
    cishu=0;
end
gama2=gama2';
area2=area2';

%% 有效反射块面积
n2=[0,0,1];
res_adjust=ans_adjust;
B3=B;
res3=res;
for i=706:-1:1  %从后往前遍历，防止数据顶替
    if abs(res_adjust(i))>=0.6
        B3(i)=[];
        res3(i,:)=[];
    end
end
u=1;
S_work=0;
for i=2:4301
    for j=1:3
        for k=1:683
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
        vec3(u,:)=n1;
        gama3(u) = dot(n1,n2)/(norm(n1)*norm(n2)); %弧度制
        S_work=S_work+Area(a,b,c)*abs(gama3(u));
        area3(u)=Area(a,b,c);
        u=u+1;
    end
    cishu=0;
end
gama3=gama3';
area3=area3';
use_old=S_jizhun/S_zong;
use_new=S_work/S_zong;


