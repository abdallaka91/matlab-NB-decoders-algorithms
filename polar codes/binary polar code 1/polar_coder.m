function u = polar_coder(msg, Q1, K)
N = length(Q1);
n = log2(N);
u = zeros(1, N);
u(Q1(N-K+1:end)) = msg;
m = 1;
for d=n-1:-1:0
    for i=1:2*m:N
        i1=i:i+m-1;
        i2=i+m:i+2*m-1;
        i3 = i:i+2*m-1;
        a = u(i1);
        b= u(i2);
        u(i3)=[mod(a+b, 2) b];

    end
    m=m*2;
end

