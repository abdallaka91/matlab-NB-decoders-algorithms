function t1=gray_code(p)

N=2^p;
t1=zeros(N,p);

for i = 1:p
    temp=[zeros(2^(i-1),1);ones(2^(i-1),1)];
    sz1=length(temp);
    while sz1<N
        temp=[temp; flip(temp)];
        sz1=length(temp);
    end
    t1(:,i)=temp;
end