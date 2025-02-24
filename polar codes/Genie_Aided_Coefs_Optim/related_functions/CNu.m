function L_theta_l = CNu(L, coefs, add_mat, div_mat)
nl1 = size(L,2);
nl=nl1/2;
q=size(L,1);
i1=1:nl;
i2=nl+1:nl1;
L_theta_l = +inf(size(L,1),nl);
L_theta_l_1=L(:,i1);
L_phi_l_1=L(:,i2);
for k=1:nl
    for i=0:q-1
        for j=0:q-1
            a_gf = add_mat(i+1,1+div_mat(j+1, 1+coefs(k)));
            a_llr = L_theta_l_1(i+1,k)+L_phi_l_1(j+1,k);
            if L_theta_l(a_gf+1,k)>a_llr
                L_theta_l(a_gf+1,k)=a_llr;
            end
        end
    end
%     [~,x1]=min(L_theta_l_1(:,k));x1=x1-1;
%     [~,x2]=min(L_phi_l_1(:,k));x2=x2-1;
%     [~,y1]=min(L_theta_l(:,k));y1=y1-1;
%     v1=add_mat(x1+1,1+div_mat(x2+1,1+coefs(k)));
%     c1=y1==v1;
end
end


