function R_Leaf = f_g_q(HD_L, Root, q, add_mat)

% v1 = Root(1,:);
% v2 = Root(2,:);
% R_Leaf = zeros(1,q);
% for i = 0 : q-1
%     i1 = add_mat(HD_L+1,i+1)+1;
%     i2 = i+1;
%     R_Leaf(i2) = v1(i1)*v2(i2);
% end
% R_Leaf = R_Leaf/sum(R_Leaf);

N = size(Root,1);
m=N/2;
R_Leaf = zeros(m,q);
ii = 1;
i1= ii:ii+m-1;
i2=ii+m:ii+2*m-1;
v1 = Root(i1,:);
v2 = Root(i2,:);
for k = 1 : m
    for i = 0 : q-1
        i1 = add_mat(HD_L(k)+1,i+1)+1;
        i2 = i+1;
        R_Leaf(k,i2) = v1(k,i1)*v2(k,i2);
    end
    R_Leaf(k,:) = R_Leaf(k,:)/sum(R_Leaf(k,:));
end


