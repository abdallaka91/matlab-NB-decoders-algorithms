clear

hamming_distance = @(a, b) arrayfun(@(x) sum(dec2bin(bitxor(a, x)) == '1'), b);

q=64;
p=log2(q);
t1d=0:q-1;
jj_max=inf;
% t2d=[0	40	19	59	38	14	53	29	15	39	28	52	41	1	58	18	30	54	13	37	56	16	43	3	17	57	2	42	55	31	36	12	60	20	47	7	26	50	9	33	51	27	32	8	21	61	6	46	34	10	49	25	4	44	23	63	45	5	62	22	11	35	24	48];
while 1
    t2d=randperm(q)-1;
    cmbs = nchoosek(1:q, 2);
    lc=size(cmbs,1);

    hmdst1=zeros(lc,1);
    for i=1:lc
        hmdst1(i)=hamming_distance(t1d(cmbs(i,1)), t1d(cmbs(i,2)));
    end

    hmdst2=zeros(lc,1);
    for i=1:lc
        hmdst2(i)=hamming_distance(t2d(cmbs(i,1)), t2d(cmbs(i,2)));
    end

    dst1=1;

    vv=[hmdst1 hmdst2];
    idx1=hmdst1==2;
    s1=sum(hmdst2(idx1));
    idx2=hmdst1==3;
    s2=sum(hmdst2(idx2));
    vv1=[hmdst1(idx1) hmdst2(idx1)];
    vv2=[hmdst1(idx2) hmdst2(idx2)];

    jj=sum(vv1(:,2)==2);
    if jj<jj_max
        jj_max=jj;
        t2d1=t2d;
    end
end