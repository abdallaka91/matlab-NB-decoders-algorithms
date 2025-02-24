function decw = polar_nb_dec(gf_coef,add_mat, mul_mat, div_mat, LLR,Ic,UIc, M,N,q,n )
% polar_nb_dec(add_mat, mul_mat, div_mat, LLR,Ic,Idxn, M,N,q,n );
L = zeros(n+1, N,q);
V = nan(n+1, N);
L(1,:,:) = LLR';
V(n+1,Ic+1)=UIc;
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

while i>0
    if stat_v{i}(j)
        if mod(j,2)==0
            i3=i-1;
            i4=n-i3+1;
            i5=Int{i3}(j/2,size(Int{i3},2)/2+1:size(Int{i3},2));
            i=i-1;
            j=j/2;
            ii1 = Int{i}(j,:);
            temp0 = sqz(V(i+1,ii1),1);

            V(i,ii1) = f_c_q(add_mat,mul_mat, temp0, gf_coef(i5, i4));
            stat_v{i}(j) = true;

        else
            i = i-1;
            j=(j+1)/2;
        end
    elseif stat_v{i+1}(2*j-1)
        ii1=Int{i}(j,:);
        ii2=Int{i+1}(2*j-1,:);
        ii3=Int{i+1}(2*j,:);
        temp0 = sqz(L(i,ii1,:),1);
        temp1 = VNu(V(i+1,ii2), temp0 , gf_coef(ii3,n-i+1), q, add_mat, mul_mat);
        i = i+1;
        j=2*j;
        ii1=Int{i}(j,:);
        L(i,ii1,:) = temp1;
        if i==n+1
            stat_v{i}(j)=true;
            if isnan(V(i,j))
                V(i,j) = HD_N_q(sqz(L(i,j,:),1));
            end
        end

    else
        ii1=Int{i}(j,:);
        ii2=Int{i+1}(2*j,:);
        temp0 = sqz(L(i,ii1,:),1);
        temp1 = CNu(temp0', gf_coef(ii2,n-i+1), add_mat, div_mat)';
        i=i+1;
        j=2*j-1;
        ii1=Int{i}(j,:);
        L(i,ii1,:) = temp1;
        if i==n+1
            stat_v{i}(j)=true;
            if isnan(V(i,j))
                V(i,j) = HD_N_q(sqz(L(i,j,:),1));
            end
            i=i-1;
            j=(j+1)/2;
        end
    end
end
decw=V(end,:);