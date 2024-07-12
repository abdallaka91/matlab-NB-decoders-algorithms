clear
pth1 = (fullfile(pwd, 'related_functions'));
addpath(pth1);
pth2 = (fullfile(pwd, 'related_variables'));
pth3 = (fullfile(pwd, 'related_variables/GF_arithm'));
pth4 = (fullfile(pwd, 'related_variables/alists'));

p=6;
q=2^p;
ex = (0:q-2)';
ex2int = [ex exponGF_to_intGF(0:q-2, p)'];
[a,b] = sort(ex2int(:,2));

int2ex = [ex2int(b,2) ex2int(b,1)];

H_matrix_mat_fl_nm = 'BeiDou_44_bb_GF64';
load([fullfile(pth4, H_matrix_mat_fl_nm) '.mat']);

[m, n] = size(h);

dv = zeros(1, n);
dc = zeros(1, m);

fst_line = [m, n, 64];
for i = 1:m
    for j = 1:n
        if h(i,j) ~= 0
            dv(j) = dv(j) + 1;
            dc(i) = dc(i) + 1;
        end
    end
end

filename = 'output.txt';
fileID = fopen(filename, 'w'); % Open the file in write mode

% Write fst_line, dc, and dv to the file
fprintf(fileID, '%d %d %d\n', fst_line);
fprintf(fileID, '%s\n', sprintf('%d ', dv));
fprintf(fileID, '%s\n', sprintf('%d ', dc));

for i = 1:m
    line = ''; % Initialize line for each iteration of i
    for j = 1:n
        if h(i,j) ~= 0
            line = [line, sprintf('%d %d ', j, int2ex(h(i,j), 2))];
        end
    end
    fprintf(fileID, '%s\n', line); % Write the line to the file
end

fclose(fileID); % Close the file
