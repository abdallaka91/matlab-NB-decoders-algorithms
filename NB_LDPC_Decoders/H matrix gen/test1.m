clear

q = 4;
p = log2(q);
dc = 3;
% H = randi([1 q-1],1, dc);
H = [ 2     1     3 ];
h = gf(H, p);
S = zeros(q, dc+1);
S(1,1)=1;
GFq = 0:q-1;

for  l0 = 1 : dc
    l = l0+1;
    for d0 = GFq
        d = d0+1;
        S(d, l) = 0;
    end
    for s0 = GFq
        s = s0+1;
        for x = GFq
            d0 = s0 + gf(H(l0),p)*gf(x,p);
            d = double(d0.x)+1;
            Sx = gf(x,p)^(sum(de2bi(x,p)));
            temp = gf(S(d,l), p) + gf(S(s,l-1), p)*Sx;
            S(d,l) = double(temp.x);


        end
    end
end




