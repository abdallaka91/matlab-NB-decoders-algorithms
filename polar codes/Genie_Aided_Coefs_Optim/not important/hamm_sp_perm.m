function [sp_frwd, sp_bck, hmdst1, hmdst2]=hamm_sp_perm(t1d, t2d)
q=length(t1d);
p=log2(q);

cmbs = nchoosek(1:q, 2);
lc=size(cmbs,1);


hmdst1=zeros(lc,1);
for i=1:lc
    hmdst1(i)=hamm_dist(t1d(cmbs(i,1)), t1d(cmbs(i,2)));
end

hmdst2=zeros(lc,1);
for i=1:lc
    hmdst2(i)=hamm_dist(t2d(cmbs(i,1)), t2d(cmbs(i,2)));
end

% vv=[hmdst1 hmdst2];
sp_frwd=zeros(p,p);
for i=1:p
    for j=1:p
        idx1=hmdst1==i;
        sp_frwd(i,j)=sum(hmdst2(idx1)==j);
    end
end

sp_bck=zeros(p,p);
for i=1:p
    for j=1:p
        idx2=hmdst2==i;
        sp_bck(i,j)=sum(hmdst1(idx2)==j);
    end
end