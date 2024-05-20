function y = f_f_q(add_mat, Root, q)

N = size(Root,1);
m=N/2;
ii = 1;
i1= ii:ii+m-1;
i2=ii+m:ii+2*m-1;
v1 = Root(i1,:);
v2 = Root(i2,:);
y = zeros(N/2, q);
for k = 1 : N/2
    v11 = v1(k, :);
    v21 = v2(k,:);
    for i = 0 : q-1
        for j= 0 : q-1

            temp = add_mat(1+i, 1+j);
            y(k, temp+1) = y(k, temp+1)+v11(i+1)*v21(j+1);
        end
    end
end
