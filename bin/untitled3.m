clear
tt1=load('C:\Users\Abdallah Abdallah\Documents\VScode\PolarizationAware_NBPC\theta_1_llr.txt');
tt2=load('C:\Users\Abdallah Abdallah\Documents\VScode\PolarizationAware_NBPC\phi_1_llr.txt');


gg0=(zeros(length(tt1)));

for i = 1 : length(tt1)
    for j = 1 : length(tt1)
        gg0(i,j)=tt1(i) + tt2(j);
    end
end


tt1=load('C:\Users\Abdallah Abdallah\Documents\VScode\PolarizationAware_NBPC\theta_1_gf.txt');
tt2=load('C:\Users\Abdallah Abdallah\Documents\VScode\PolarizationAware_NBPC\phi_1_gf.txt');


gg=(zeros(length(tt1)));

for i = 1 : length(tt1)
    for j = 1 : length(tt1)
        temp = gf(tt1(i), 6) + gf(tt2(j), 6);
        gg(i,j)=double(temp.x);
    end
end

gg0a = gg0(:);
[a,b] = sort(gg0a);
gga=gg(:);
c=gga(b);
c(1:6)


