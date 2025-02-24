function R_Leaf = VNu(HD_L, Root,coefs, q, add_mat, mul_mat)

N = size(Root,1);
m=N/2;
R_Leaf = zeros(m,q);
i1= 1:m;
i2=m+1:N;
v1 = Root(i1,:);
v2 = Root(i2,:);
for k = 1 : m
    for i = 0 : q-1
        a=add_mat(HD_L(k)+1,i+1);
        b=mul_mat(i+1,coefs(k)+1);
        R_Leaf(k,i+1) = v1(k,a+1)+v2(k,b+1);
    end
    R_Leaf(k,:) = R_Leaf(k,:)-min(R_Leaf(k,:));
end
