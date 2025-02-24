function [G, G_arr, Ginv, Ginv_arr] = coefs_2_G_matrix(N, coefs, add_mat, mul_mat, div_mat)
G = ones(N,N);
G(:,1)=1;
n = round(log2(N));
h=nan(N, n);
t=0:N/2-1;
for l = 1 : n
    f=2^(l-1)+2*t-mod(t, 2^(l-1))+1;
    h(f,l) = coefs(:,l);
end

for i = 2 :2: N
    G(i,2)=coefs(i/2);
end

for l = 1 :1: n-1
    j=2^l;
    j1 = j+1;
    j2=2^(l+1);
    Nj=N/j;
    for p = 1 : Nj
        i1=(p-1)*j+1;
        i2=p*j;
        if mod(p,2)~=0
            G(i1:i2,j1:j2)=0;
        else
            a = G(i1:i2,1:j);
            for z=1:size(a,2)
                i3=i1+z-1;
                for z1 = 1 : size(a,1)
                    a(z1,z)=mul_mat(a(z1,z)+1,h(i3,l+1)+1);
                end
            end
            G(i1:i2,j1:j2)=a;
        end
    end
end
[row, col] = find(G);
G_arr=zeros(3*length(row),1);
for k = 1:length(row)
    G_arr(3*(k-1)+1) =row(k);
    G_arr(3*(k-1)+1+1) = col(k);
    G_arr(3*(k-1)+1+2) = G(row(k), col(k));
end
Ginv = matrix_inversion_GF(G, mul_mat, add_mat, div_mat);
[row, col] = find(Ginv);
Ginv_arr=zeros(3*length(row),1);
for k = 1:length(row)
    Ginv_arr(3*(k-1)+1) =row(k);
    Ginv_arr(3*(k-1)+1+1) = col(k);
    Ginv_arr(3*(k-1)+1+2) = Ginv(row(k), col(k));
end

end

