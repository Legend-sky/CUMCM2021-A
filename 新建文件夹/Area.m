function S=Area(A,B,C)  %利用余弦定理求解三角形面积
if length(B)==3  %输入三点是三维空间坐标
    AB=B-A;
    BC=C-B;
end
    Z=cross(AB,BC);  %叉乘
    S=1/2*norm(Z);   %面积（余弦定理）
end