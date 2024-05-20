function [cword, msg] = gen_rand_code(N,K, Q1)
n = log2(N);
msg =randi([0 1],1,K);
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
        cword(i3)=[mod(a+b, 2) b]; %combining

    end
    m=m*2;
end
