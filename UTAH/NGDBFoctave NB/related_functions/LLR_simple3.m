
function LLR = LLR_simple3(I,LLRfact , unreliable_sat, q,N, alph_bin, LLR)
for n = 1 : N
    mx1 = -inf;
    for i = 1 : q
        LLR(n,i) =  -round(LLRfact*sum(alph_bin(i,:).*I(n,:))); 
        if LLR(n,i)>mx1
            mx1 = LLR(n,i);
        end
    end
    LLR(n,:) = LLR(n,:) - mx1;
    LLR(n,LLR(n,:)<unreliable_sat) = unreliable_sat;
end

