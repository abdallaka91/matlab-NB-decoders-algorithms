clear
mat = [1 0 6 7 11 18 16 44 21 53  
2 0 7 7 12 18 17 44 22 53  
3 7 8 18 13 0 18 44 23 53  
4 7 9 18 14 44 19 0 24 53  
5 7 10 18 15 44 20 53 25 0  
1 7 7 0 13 18 19 44 25 53  
2 7 8 44 14 53 20 18 21 0  
3 0 9 7 15 18 16 53 22 44  
4 0 10 7 11 53 17 18 23 44  
5 18 6 44 12 53 18 7 24 0 ];

loc = mat(:,1:2:end);

values = mat(:,2:2:end);

% values here is the power of alpha
%%
p = 2; % Or any prime number
m = 6; % Or any positive integer
% field = gftuple(values(:),m,p);
% field1 = bi2de(field);
% values1 = reshape(field1,size(values,1),size(values,2));
y = alpha2dec( m, [1 1 0 0 0 0 1]);

%%
M = 10;
N = 25;
idx = loc;

vl = values;

h = (zeros(M,N));
for i = 1 : M
    lst_idx = idx(i,:);
    lst_vl = vl(i,:);

    for j = 1 : length(lst_idx)
        tmp = y(lst_vl(j)+2);
        h(i, lst_idx(j)) = tmp;
    end

end

h1 = 0*full(h);
figure;
hold on
for i = 1 : M
    jj=find(h(i,:));
    h1(i,jj)=1;
    plot(jj, i*ones(size(jj)), 'b.')
end
ff=sum(h1,1);


pth4 = (fullfile(fileparts(pwd), 'related_variables/alists/'));

% fl1 = ['generated_' num2str(M) 'x' num2str(N) '_GF' num2str(q) '.mat'];

fl1 = 'LabStics_GF64_K90_Dc5_1.mat';
fll_nm = fullfile(pth4, fl1);
save(fll_nm,'h')