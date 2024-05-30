clear
a = [4 27   7 1   10 36   16 49   
4 60   5 45   9 32   15 23   
2 6   6 21   7 56   14 47   
1 32   8 23   11 60   16 45   
3 52   9 61   11 26   14 11   
1 39   6 2   12 24   15 11   
2 36   5 45   8 10   13 58   
3 62   10 12   12 53   13 27];
%%
M = size(a,1);
idx = a(:,1:2:end);
N  = max(max(idx));

vl = a(:,2:2:end);

h = (zeros(M,N));
for i = 1 : M
    lst_idx = idx(i,:);
    lst_vl = vl(i,:);
    e1 = lst_idx == 0;
    lst_idx(e1) = [];
    lst_vl(e1) = [];
    h(i, lst_idx) = lst_vl;

end

pth4 = (fullfile(fileparts(pwd), 'related_variables\alists\'));

% fl1 = ['generated_' num2str(M) 'x' num2str(N) '_GF' num2str(q) '.mat'];

fl1 = 'N96_K48_GF64.mat';
fll_nm = fullfile(pth4, fl1);
save(fll_nm,'h')