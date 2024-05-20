function [cword, msg] = gen_rand_nb_code(N,K, Q1, q, add_mat, mul_mat, msg)


n = log2(N);
if nargin<7
    msg =randi([0 q-1],1,K);
elseif length(msg)~=K
    msg =randi([0 q-1],1,K);
end
cword = zeros(1, N);
cword(Q1(N-K+1:end)) = msg;
m = 1;
for d=n-1:-1:0
    for i=1:2*m:N
        i1=i:i+m-1;
        i2=i+m:i+2*m-1;
        i3 = i:i+2*m-1;
        a = cword(i1);
        b= cword(i2);
        %         cword(i3)=[mod(a+b, 2) b]; %combining
        c = a;
        for j = 1 : length(i1)
            c(j)=sum_arr_gf_dec([a(j) b(j)], add_mat);
        end
        cword(i3) = [c b];


    end
    m=m*2;
end
