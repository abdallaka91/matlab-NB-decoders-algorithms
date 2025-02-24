clear

% rng(0)
filePath = 'C:\Users\Abdallah Abdallah\Documents\Personals\matlab-NB-decoders-algorithms\polar codes\december 2024\bpsk_ccsk_gf64\biawgn_ccsk[6]_gfdim6_n10.plr';
[p, data] = read_nb_polar_files_parameter(filePath);
% best channels to worst channels indexes
q=2^p;
mainPath = pwd;
% mainPath = fileparts(projectPath);
related_variables_pth = fullfile(mainPath, 'related_variables');
pth3 = (fullfile(related_variables_pth, 'GF_arithm'));
fl_nm = ['arith_' num2str(q) '.mat'];
load(fullfile(pth3, fl_nm));

N=length(data(1).Pe);

for i2 = 1 : length(data)
    coefs = data(i2).gf_coef;
    %
    g1=[1 0; 1 1];
    g=g1;
    n=log2(N);
    for i= 1 : n-1
        g1=kron(g1,g);
    end
    G=zeros(N,N);
    for i = 1 : N
        i1 = find(g1(:,i));
        N1=length(i1);
        succss=false;
        Tries=0;
        while ~succss && Tries<7
            Tries=Tries+1;
            a=zeros(N1);
            b=zeros(N1,1);
            for j=1:N1
                c1=true;
                while c1
                    x=randi([0 q-1],1,N);
                    y=encode(x, coefs,add_mat, mul_mat);
                    a(j,:)=x(i1);
                    b(j)=y(i);
                    if i==N
                        if b(j)~=0
                            c1=false;
                        else
                            fprintf('Null!\n')
                        end
                    else
                        c1=false;
                    end
                end
            end
            try
                a1=matrix_inversion_GF(a, mul_mat, add_mat, div_mat);
                succss=true;
            catch me
                fprintf('Attempt %d failed (@ i = %d)\n', Tries, i);
            end
        end

        c=gf_mat_mul(a1,b,add_mat, mul_mat);
        G(i1,i)=c;
    end
    %%

    [row, col] = find(G);
    G_arr=zeros(3*length(row),1);
    for k = 1:length(row)
        G_arr(3*(k-1)+1) =row(k);
        G_arr(3*(k-1)+1+1) = col(k);
        G_arr(3*(k-1)+1+2) = G(row(k), col(k));
    end
    Ginv = matrix_inversion_GF(G, mul_mat, add_mat, div_mat);
    [row, col] = find(Ginv);
    Ginv_arr=zeros(3*length(row),1);
    for k = 1:length(row)
        Ginv_arr(3*(k-1)+1) =row(k);
        Ginv_arr(3*(k-1)+1+1) = col(k);
        Ginv_arr(3*(k-1)+1+2) = Ginv(row(k), col(k));
    end

    data(i2).G=G;
    data(i2).G_arr=G_arr;
    data(i2).Ginv=Ginv;
    data(i2).Ginv_arr=Ginv_arr;

    x=randi([0 q-1],1,N);
    y=encode(x, coefs,add_mat, mul_mat);
    y1 = GF_pol_enc_dec_mat_mul(x,G_arr, mul_mat, add_mat);
    x1 = GF_pol_enc_dec_mat_mul(y1,Ginv_arr, mul_mat, add_mat);
    p1=prod(y==y1);
    p2=prod(x==x1);
    pp=p1*p2;
    if ~pp==1
        fprintf('Errooooooooor!\n')
    end
end
[~, fileName, ~] = fileparts(filePath);

folder_of_mat_file = fullfile(pwd, 'coefs_G_Ginv');

save(fullfile(folder_of_mat_file, [fileName '.mat']), 'data')



