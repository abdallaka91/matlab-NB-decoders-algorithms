clear
pth1 = (fullfile(pwd, 'related_functions'));
q=4;

add_mat = GF_arithm_matrix(q, 'add');
mul_mat = GF_arithm_matrix(q, 'mul');
div_mat = GF_arithm_matrix(q, 'div');

addpath(pth1);
h=[1 2 1 0 0 0; ...
    1 0 0 3 0 1; ...
    0 0 3 2 2 0; ...
    0 1 0 0 3 2];


[G,~] = Generator_matrix_G_from_full_rank_H(h, add_mat, mul_mat, div_mat);
info_seq =[ 1 3];
code_seq = gf_mat_mul(info_seq,G, add_mat, mul_mat);
valid_symdrom = gf_mat_mul(code_seq,h', add_mat, mul_mat);

h0=h;
for i = 1 : 4
    str_cn_vn0{i}=find(h0(i,:));
end

h(1,:) = gf_mat_mul(div_mat(1+1, h(1,1)+1),h(1,:),add_mat,mul_mat);
h(2,:) = gf_mat_mul(div_mat(1+1, h(2,1)+1),h(2,:),add_mat,mul_mat);
h(2,:)= gf_mat_add(h(2,:), h(1,:), add_mat);

h(2,:) = gf_mat_mul(div_mat(1+1, h(2,2)+1),h(2,:),add_mat,mul_mat);
h(3,:) = gf_mat_mul(div_mat(1+1, h(3,3)+1),h(3,:),add_mat,mul_mat);

h(4,:)= gf_mat_add(h(4,:), h(2,:), add_mat);
h(4,:) = gf_mat_mul(div_mat(1+1, h(4,3)+1),h(4,:),add_mat,mul_mat);
h(4,:)= gf_mat_add(h(4,:), h(3,:), add_mat);
h(4,:) = gf_mat_mul(div_mat(1+1, h(4,5)+1),h(4,:),add_mat,mul_mat);

h(2,:) = gf_mat_mul(div_mat(1+1, h(2,6)+1),h(2,:),add_mat,mul_mat);
h(2,:)= gf_mat_add(h(2,:), h(4,:), add_mat);

h(2,:)= gf_mat_add(h(2,:), h(3,:), add_mat);
h(1,:)= gf_mat_add(h(1,:), h(2,:), add_mat);

for i = 1 : 4
    str_cn_vn{i}=find(h(i,:));
end

v=zeros(1,6);
vv=zeros(16,6);
vvv = zeros(16,4);
k=1;
for i = 0 : 3
    for j = 0 : 3
        v([3 5]) = [i j];
        v(1) = add_mat(v(3)+1, mul_mat(2+1, v(5)+1)+1);
        v(2)=v(5);
        v(4) = div_mat(add_mat(v(3)+1, mul_mat(v(5)+1, 3+1)+1)+1, 3+1);
        v(6) = v(5);
        vv(k,:)=v;
        vvv(k,:)=v([1 2 4 6]);
        k=k+1;
        y0 = decod_prod(v,h0,str_cn_vn0, mul_mat, add_mat)
    end

end
