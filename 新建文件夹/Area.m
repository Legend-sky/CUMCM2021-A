function S=Area(A,B,C)  %�������Ҷ���������������
if length(B)==3  %������������ά�ռ�����
    AB=B-A;
    BC=C-B;
end
    Z=cross(AB,BC);  %���
    S=1/2*norm(Z);   %��������Ҷ���
end