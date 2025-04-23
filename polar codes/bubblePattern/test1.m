clear
aa=load("bubbles_N64_GF64_SNR-7.50_18x18_Cs_mat.txt");

Nobs = 10000;
N=64;
n=log2(N);
q=64;
nL = size(aa,2);
nH = size(aa,1)/sum(2.^(0:n-1));
Pt = 0.12;

Cll = cell(n,1);
Cll_mat = cell(1,n);
Cll_norm_cll = cell(1,n);
Cll_norm_mat = cell(1,n);

Cll_bin = cell(1,n);
k=0;
for i = 0:n-1
    Cll{i+1} = cell(2^i,1);
    Cll_mat{i+1} = zeros(nH, nL, 2^i);
    Cll_norm_mat{i+1} = zeros(nH, nL, 2^i);
    Cll_norm_cll{i+1} = cell(1,2^i);
    Cll_bin{i+1} = zeros(nH, nL, 2^i);
    for j=0:2^i-1
        Cll{i+1}{j+1}= aa(k*nH+1:(k+1)*nH,:);
        Cll_mat{i+1}(:,:,j+1) = Cll{i+1}{j+1};
        Cll_norm_mat{i+1}(:,:,j+1) = Cll_mat{i+1}(:,:,j+1);%/(2^(n-(i+1)))/Nobs;
        Cll_norm_cll{i+1}{j+1} = Cll_norm_mat{i+1}(:,:,j+1);
        Cll_bin{i+1}(:,:,j+1) = Cll{i+1}{j+1}==0;
        k=k+1;
    end
end


indic_mat_idx = cell(1,n);
indic_cell = Cll_norm_cll;
indic_mat = Cll_norm_mat;
for i = 0:n-1
    indic_mat_idx{i+1} = cell(2^i,1);
    for j=0:2^i-1
        indic_cell{i+1}{j+1}=indic_cell{i+1}{j+1}>Pt;
        indic_mat{i+1}(:,:,j+1)=squeeze(indic_mat{i+1}(:,:,j+1))>Pt;
        indic_mat_idx{i+1}{j+1} = [];
        temp_mat = Cll_norm_mat{i+1}(:,:,j+1);
        for i1 = 1 : nH
            for j1 = 1 : nL
                if(temp_mat(i1,j1))>Pt
                    indic_mat_idx{i+1}{j+1}=[indic_mat_idx{i+1}{j+1}; [i1-1 j1-1]];
                end
            end
        end
    end
end
plot2D_bubble(Cll_norm_mat, nH, nL, n);
plot2D_bubble(indic_mat, nH, nL, n);
