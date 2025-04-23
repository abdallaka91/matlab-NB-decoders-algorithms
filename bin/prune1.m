clear
N=64;
n=log2(N);

cn1 = {
 14; 
[18 4]; 
[0 10 5 1]; 
[0 0 11 2 8 1 1 1]; 
[0 0 0 7 0 4 3 1 0 1 1 1 1 1 1 1]; 
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]}; 

for l=n-1:-1:2
    for s = 1:2:length(cn1{l})
        ln1 = max(cn1{l}(s), cn1{l-1}((s+1)/2));
        cn1{l-1}((s+1)/2)=ln1;
    end
end

for i = 1 : n-1
    disp(cn1{i})
end