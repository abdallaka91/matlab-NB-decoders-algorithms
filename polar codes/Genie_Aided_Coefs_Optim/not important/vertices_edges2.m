clear
hamming_distance = @(a, b) arrayfun(@(x) sum(dec2bin(bitxor(a, x)) == '1'), b);

p=6;
np=2;
perm_i = [np p-np];
q=2^p;
dsts=zeros(p,1);
nb_edges=dsts;
words1=(0:q-1)';

for i = 1 :p
    dsts(i)=nchoosek(p,i);
    nb_edges(i) =dsts(i)*q/2;
end
lst1 = nan(q,2);
lst1_bin = (repmat("*",q,1));

pp1=0;

for dd=perm_i
    cnt_filledp=inf;
    pp1=pp1+1;
    edges_nodes = zeros(nb_edges(perm_i(1)),2);
    nodes_edges = zeros(q,dsts(perm_i(1)));
    nodes_edges_bin = string(nodes_edges);

    for i=1:q
        a=string(dec2bin(i-1, p));
        [~, lst2b] = hmm3([a a],p, [dd dd]);
        for j=1:dsts(perm_i(1))
            nodes_edges_bin(i,j)=string(lst2b(j,:));
            nodes_edges(i,j)=bin2dec((lst2b(j,:)));
        end
    end

    vert_labl = (0:q-1)';
    vert_val_bin = (repmat("*",q,1));

    vert_val_bin(1)=string(dec2bin(0,p));

    cnt_filled=1;

    for i=1:2%dsts(dd)
        vert_val_bin(i+1,1)=nodes_edges_bin(1,i);
        cnt_filled=cnt_filled+1;
    end

    cnt_filled1=0;
    while cnt_filled<q
        cmb_vert = (nchoosek(1:cnt_filled,2));
        cmb_cnt = size(cmb_vert,1);
        for j=1:cmb_cnt


            [~, lst2b] = hmm3(vert_val_bin(cmb_vert(j, :)),p, [dd dd]);
            lst2b=string(lst2b);
            aaa=1;

            if size(lst2b,1)>1
                for i=1:size(lst2b,1)
                    k=0;
                    ext=0;
                    while k<=cnt_filled
                        k=k+1;
                        if(strcmp(lst2b(i), vert_val_bin(k)))
                            ext=1;
                            break;
                        end
                    end
                    if ext==0
                        cnt_filled=cnt_filled+1;
                        vert_val_bin(cnt_filled)=lst2b(i);
                    end
                end
            end
        end
        if cnt_filled1==cnt_filled
            break
        else
            cnt_filled1=cnt_filled;
        end
    end
    if cnt_filled<cnt_filledp
        cnt_filledp=cnt_filled;
    end
    for i =1:cnt_filled
        lst1(i,pp1)=bin2dec(vert_val_bin(i));
        lst1_bin(i,pp1)=(vert_val_bin(i));
    end


    if cnt_filled1~=q
        q1=cnt_filled1;
        for i=1:2%dsts(dd)
            vert_val_bin(cnt_filled+1,1)=nodes_edges_bin(2,i);
            cnt_filled=cnt_filled+1;
        end
    end

    cnt_filledp=inf;
    while cnt_filled<q
        if cnt_filled1==cnt_filled
            break
        else
            cnt_filled1=cnt_filled;
        end
        cmb_vert = (nchoosek(1:cnt_filled,2));
        cmb_cnt = size(cmb_vert,1);
        for j=1:cmb_cnt


            [~, lst2b] = hmm3(vert_val_bin(cmb_vert(j, :)),p, [dd dd]);
            lst2b=string(lst2b);
            aaa=1;

            if size(lst2b,1)>1
                for i=1:size(lst2b,1)
                    k=0;
                    ext=0;
                    while k<=cnt_filled
                        k=k+1;
                        if(strcmp(lst2b(i), vert_val_bin(k)))
                            ext=1;
                            break;
                        end
                    end
                    if ext==0
                        cnt_filled=cnt_filled+1;
                        vert_val_bin(cnt_filled)=lst2b(i);
                    end
                end
            end
        end
    end
    if cnt_filled<cnt_filledp
        cnt_filledp=cnt_filled;
    end
    for i =1:cnt_filled
        lst1(i,pp1)=bin2dec(vert_val_bin(i));
        lst1_bin(i,pp1)=(vert_val_bin(i));
    end

end


%%
lst1=lst1(1:cnt_filledp,:);
lst1_bin=lst1_bin(1:cnt_filledp,:);
q=cnt_filledp;
% [~,idxs]=sort(lst1(:,1));
% lst1=lst1(idxs,:);
% lst1_bin=lst1_bin(idxs,:);

t1d=lst1(:,1);
t2d=lst1(:,2);

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

idx1=hmdst1==(perm_i(1));
s1=sum(hmdst2(idx1));
idx2=hmdst1==perm_i(2);
s2=sum(hmdst2(idx2));
vv1=[hmdst1(idx1) hmdst2(idx1)];
vv2=[hmdst1(idx2) hmdst2(idx2)];