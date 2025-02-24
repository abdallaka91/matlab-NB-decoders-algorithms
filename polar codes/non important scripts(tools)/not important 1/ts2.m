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

%
g1=[1 0; 1 1];
g=g1;
n=log2(N);
for i= 1 : n-1
    g1=kron(g1,g);
end
a2=cell(N,1);
b2=a2;
f2=a2;
l2=zeros(N,1);
p2=l2;
cnt2=ones(N,1);

for i = 1 : N
    f2{i}=find(g1(:,i));
    l2(i)=length(f2{i});
    a2{i}=zeros(l2(i));
    b2{i}=zeros(l2(i),1);
end
i2 = 1;
coefs = data(i2).gf_coef;
G=zeros(N,N);
dne=false;
Tries=0;
while dne==false && Tries<5
    Tries=Tries+1;
    cnt2=ones(N,1);
    MM=max(l2(p2==0));
    for i = 1 : MM
        x=randi([0 q-1],1,N);
        y=encode(x, coefs,add_mat, mul_mat);
        for ii=1:N
            if p2(ii)==0
                if cnt2(ii)<=l2(ii)
                    a2{ii}(cnt2(ii),:)=x(f2{ii});
                    b2{ii}(cnt2(ii))=y(ii);
                    cnt2(ii)=cnt2(ii)+1;
                end
            end
        end
    end
    for i=1:N
        if p2(i)==0
            a=a2{i};
            b=b2{i};
            try
                a1=matrix_inversion_GF(a, mul_mat, add_mat, div_mat);
                c=gf_mat_mul(a1,b,add_mat, mul_mat);
                G(f2{i},i)=c;
                p2(i)=true;
            catch me
                fprintf('failed(Attempt = %d, i = %d)\n',Tries, i);
            end
        end
    end
    if prod(p2)==1
        dne=true;
    end
end

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

