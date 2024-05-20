function prb_q = prb_q(y_nse, sigm)
p = size(y_nse,2);
q=2^p;
N = size(y_nse,1);
I = 2*y_nse/sigm^2;
alph = 0:q-1;
alph_bin = fliplr(dec2bin(alph, p) - 48);
prb1 = y_nse;
prb0 = y_nse;

prb_q = zeros(N,q);
for i = 1 : N
    for j = 1 : p
        prb1(i,j) = 1/(1+exp(I(i,j)));
    end
end
prb0 = 1 - prb1;
for i = 1 : N
    for j = 1 : q
        temp_prb = 1;
        for k = 1 : p
            if alph_bin(j,k)==1
                temp_prb = temp_prb*prb1(i,k);
            else 
                temp_prb = temp_prb*prb0(i,k);
            end
        end
        prb_q(i,j) = temp_prb;
    end
end
end
                






