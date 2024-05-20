function decw = polar_nb_dec(add_mat, prb,frozen_idx,M,N,q,n )

L = zeros(n+1, N,q);
V = zeros(n+1, N);
L(1,:,:) = prb;

stat_v = cell(n+1,1);
for i = 1 : n+1
    i0 = i-1;
    stat_v{i,1} = false(1,  2^(i0));
end


I = cell(n,1);
Int = cell(n,1);
for i = 1 : n+1
    mi = 2^(i-1);
    li = N/mi;
    for j = 0:mi-1
        I{i}(j+1) =j*li+1;
        Int{i}(j+1,:) = I{i}(j+1):I{i}(j+1)+li-1;
        
    end
end


i = 1;
j = 1;
froz_j1 = 1;
decw = nan(1,N);
while i>0
    if stat_v{i}(j)
        if mod(j,2)==0
            i=i-1;
            j=j/2;
            ii1 = Int{i}(j,:);
            temp0 = sqz(V(i+1,ii1),1);
            V(i,ii1) = f_c_q(add_mat, temp0);
            stat_v{i}(j) = true;

        else
            i = i-1;
            j=(j+1)/2;
        end
    elseif stat_v{i+1}(2*j-1)
        ii1=Int{i}(j,:);
        temp0 = sqz(L(i,ii1,:),1);
        ii1=Int{i+1}(2*j-1,:);
        temp1 = f_g_q(V(i+1,ii1), temp0 , q, add_mat);
        i = i+1;
        j=2*j;
        ii1=Int{i}(j,:);
        L(i,ii1,:) = temp1;
        if i==n+1
            stat_v{i}(j)=true;
            V(i,j) = HD_N_q(sqz(L(i,j,:),1));
            if froz_j1<=M
                if j==frozen_idx(froz_j1)
                    V(i,j) = 0;
                    froz_j1 = froz_j1+1;
                end
            end
            decw(1:j) = V(i,1:j);
        end
    
    else
        
        temp0 = sqz(L(i,Int{i}(j,:),:),1);
        temp1 = f_f_q(add_mat, temp0, q);
        i=i+1;
        j=2*j-1;
        ii1=Int{i}(j,:);
        L(i,ii1,:) = temp1;
        if i==n+1
            stat_v{i}(j)=true;
            V(i,j) = HD_N_q(sqz(L(i,j,:),1));
            if froz_j1<=M
                if j==frozen_idx(froz_j1)
                    V(i,j) = 0;
                    froz_j1 = froz_j1+1;
                end
            end
            decw(1:j) = V(i,1:j);
            i=i-1;
            j=(j+1)/2;
        end
    end
end