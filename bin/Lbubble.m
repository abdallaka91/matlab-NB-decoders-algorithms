function [E, E_gf] = Lbubble(u,v,u_gf, v_gf, add_mat, nL, nH, nopM, nb, nc, zC)


if u(zC)>v(zC)
    a = u;a_gf = u_gf;
    b = v; b_gf = v_gf;
else
    a = v;a_gf = v_gf;
    b = u; b_gf = u_gf;
end

T = nan(nH, nL);
T_gf = nan(nH,nL);

cnt=0;
E=[];
E_gf=[];
for i = 1 : nb
   T_gf(i,1) = add_mat(a_gf(i)+1, b_gf(1)+1);
end
T(1:nb,1) = a(1:nb);
srtr = T(1:nb,1);
srtrg = T_gf(1:nb,1);
sij = [1:nb; ones(1, nb)];


Ti = false(size(T));
Ti(1:nb,1)=true;
nop = 0;
while nop<nopM
    nop=nop+1;
    [m,n] = min(srtr);
    i=sij(1,n);
    j=sij(2,n);
    if sum(E_gf==srtrg(n))==0
        E=[E m];
        E_gf=[E_gf srtrg(n)];
        cnt = cnt+1;
    end
    if i==nH ||j==nL || cnt==nc
        break
    end
    if i==1
        H = 1; Hb=~H;
    end
    if i==nb
        H = 0; Hb=~H;
    end
    if ~Ti(i+Hb, j+H)
        i1 = i+Hb;
        j1 = j+H;
    else
        i1 = i+H;
        j1 = j+Hb;
    end
    Ti(i1, j1) = true;
    T(i1,j1) = a(i1)+b(j1);
    T_gf(i1,j1) = add_mat(a_gf(i1)+1, b_gf(j1)+1);
    srtr(n) = T(i1,j1);
    srtrg(n) = T_gf(i1,j1);
    sij(1,n) = i1;
    sij(2,n) = j1;
    if  nop==nopM
        break
    end
end

