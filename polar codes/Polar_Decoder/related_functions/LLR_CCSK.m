function L  = LLR_CCSK(y, etaq, q, N, sigma2)
L = zeros(size(y));
HD1 = y<0;
for n = 1 : N
    for j = 1 : q
        temp = (2*y(:,n)/(sigma2^2)).*(etaq(:,j)-HD1(:,n));
        L(j,n)=sum(temp);
    end
%     L(:,n)=L(:,n)-min(L(:,n));
end
end

