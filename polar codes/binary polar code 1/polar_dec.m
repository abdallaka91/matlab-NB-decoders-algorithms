function decw = polar_dec( rn, L, stat_v, V,frozen_idx,M,n )

for i = 1 : n+1
    i0=i-1;
    stat_v{i,1} = false(1,  2^(i0));
end
L{1} = rn;

i = 1;
j = 1;
froz_j1 = 1;
decw = nan(1,length(L{1}));
while i>0
    if stat_v{i}(j)
        if mod(j,2)==0
            i=i-1;
            j=j/2;
            temp0 = V{i+1}(2*j-1,:);
            temp1 = V{i+1}(2*j,:);
            temp = function_comb(temp0, temp1);
            V{i}(j,:) = temp;
            stat_v{i}(j) = true;

        else
            i = i-1;
            j=(j+1)/2;
        end
    elseif stat_v{i+1}(2*j-1)
        i = i+1;
        j=2*j;
        temp0 = L{i-1}(j/2,:);
        temp1 = V{i}(j-1,:);
        temp = function_g(temp0, temp1);
        L{i}(j,:) = temp;
        if i==n+1
            stat_v{i}(j)=true;
            Li = L{i}(j,:);
            V{i}(j,1) = double(Li<0);
            if froz_j1<=M
                if j==frozen_idx(froz_j1)
                    V{i}(j,1) = 0;
                    froz_j1 = froz_j1+1;
                end
            end
            decw(1:j) = V{i}(1:j);
        end
    else
        i=i+1;
        j=2*j-1;
        temp0 = L{i-1}((j+1)/2,:);
        temp = function_f(temp0);
        L{i}(j,:) = temp;
        if i==n+1
            stat_v{i}(j)=true;
            Li = L{i}(j,:);
            V{i}(j,1) = double(Li<0);
            if froz_j1<=M
                if j==frozen_idx(froz_j1)
                    V{i}(j,1) = 0;
                    froz_j1 = froz_j1+1;
                end
            end
            decw(1:j) = V{i}(1:j);
            i=i-1;
            j=(j+1)/2;
        end
    end
end