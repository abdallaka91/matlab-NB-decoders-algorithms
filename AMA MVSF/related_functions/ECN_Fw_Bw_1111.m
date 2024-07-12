function [c2v_llr, c2v_gfp] = ECN_Fw_Bw_1111(v2c_gfp,v2c_llr, dc,nm, add_mat)
c2v_llr=[];
c2v_gfp=[];
u=1;
% clear
% q=16;
% nm = q;
% dc = 6;
%     add_mat = GF_arithm_matrix(q, 'add');
%     mul_mat = GF_arithm_matrix(q, 'mul');
%     div_mat = GF_arithm_matrix(q, 'div');
% v2c_llr = randn(dc,nm);
% [v2c_llr, b]=sort(v2c_llr,2);
% v2c_gfp = b-1;
% nm = size(v2c_gfp,2);
if u
    v2c_llr=bsxfun(@minus, v2c_llr, min(v2c_llr')');
    cmb1 = 1:nm;
    cmb = cmb1;
    for i = 1 : dc-1
        cmb = CombVec(cmb1, cmb);
    end
    cmb=cmb';
    dlt_gf = cmb;
    dlt_llr = dlt_gf;
    L = size(cmb,1);
    snd_llr = zeros(L,1);
    snd_gf = snd_llr;

    for i = 1 : L
        for j = 1 : dc
            snd_llr(i) = snd_llr(i)+v2c_llr(j, cmb(i, j));
            snd_gf(i) = add_mat(snd_gf(i)+1, v2c_gfp(j, cmb(i, j))+1);
            dlt_gf(i, j) = v2c_gfp(j, cmb(i, j));
            dlt_llr(i, j) = v2c_llr(j, cmb(i, j));
        end
    end
    %%
    [snd_llr_sort, b] = sort(snd_llr);
    snd_gf_sort = snd_gf(b);
    clear snd_llr snd_gf cmb
    dlt_gf_sort = dlt_gf(b, :);
    clear dlt_gf
    dlt_llr_sort = dlt_llr(b,:);
    clear dlt_llr
    c2v_gfp = nan(dc, nm);
    c2v_llr = c2v_gfp;
    kk = zeros(1,dc);
    for i = 1 : L
        for j = 1 : dc
            v1 = add_mat(snd_gf_sort(i)+1, dlt_gf_sort(i,j)+1);
            l1 = snd_llr_sort(i)-dlt_llr_sort(i, j);
            if sum(c2v_gfp(j,:)==v1)==0 && kk(j)<nm
                kk(j)=kk(j)+1;
                c2v_gfp(j,kk(j))=v1;
                c2v_llr(j,kk(j))=l1;
            end
            if prod(kk==nm)==1
                break
            end
        end
    end
    clear b snd_llr_sort dlt_gf_sort dlt_llr_sort snd_gf_sort
end
