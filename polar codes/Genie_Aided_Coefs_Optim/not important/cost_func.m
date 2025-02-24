function [cnt, sp_frwd]= cost_func(t2d, wght)
t2d=t2d-1;
q=length(t2d);
t1d=0:q-1;
cmbs = nchoosek(1:q, 2);
lc=size(cmbs,1);

hmdst1=zeros(lc,1);
for i=1:lc
    hmdst1(i)=hamm_dist(t1d(cmbs(i,1)), t1d(cmbs(i,2)));
end

[sp_frwd, sp_bck, hmdst1, hmdst2]=hamm_sp_perm(t1d, t2d);

arr1=[sp_frwd(1,1) sp_frwd(1,2) sp_frwd(2,1) sp_frwd(2,2) sp_frwd(1,3) sp_frwd(3,1)];
arr2=arr1.*wght;
cnt = sum(arr2);

    
end


