clear
Nf = 1e5;
sgm = 0.5;
nse = randn(2,Nf)*sgm;
frms = randi([0 1],2,Nf);
frms(1,:)=0;
cd_frms = frms;
for i = 1:Nf
    m = frms(:,i);
    CD = [xor(m(1), m(2)) m(2)];
    cd_frms(:,i) =CD;
end
md_cd_frms = -2*cd_frms+1;
y = md_cd_frms+nse;

HD = double(y<0);

df_hd = double(HD~=cd_frms);
ndf_hd=sum(df_hd,2);

I = 2*y/(sgm^2);
prb1 = 1./(1+exp(I));
prb0 = exp(I)./(1+exp(I));

p0u0 = prb0(1,:).*prb0(2,:) + prb1(1,:).*prb1(2,:);
p0u1 = prb0(1,:).*prb1(2,:) + prb1(1,:).*prb0(2,:);
u0 = 1*(p0u0<p0u1);

p1u0 = p0u0*0;
p1u1 = p0u0*0;
for i = 1 : Nf

%     if u0(i)<0.5
        p1u0(i) = prb0(1,i)*prb0(2,i);
        p1u1(i) = prb1(1,i)*prb1(2,i);
%     else
%         p1u0(i) = prb1(1,i)*prb0(2,i);
%         p1u1(i) = prb0(1,i)*prb1(2,i);
%     end
end

u1 = p1u0<p1u1;


ndf1 = sum(u0~=frms(1,:))
ndf2 = sum(u1~=frms(2,:))



