function e=cond_entropy(p)
p(p<eps)=eps;
l1=log2(1./p);
p1=p.*l1;
e=sum(p1,1);
end

