function save_idx_val_PCM(file_name, N,M,q, idx_val, dc, dv)

clear conf_detail
conf_detail.a1frstline = sprintf("%d %d %d", N, M, q);
conf_detail.a2scndline = sprintf("%d ",repmat(dv,1,N));
conf_detail.a3scndline = sprintf("%d ",repmat(dc,1,M));
msg='';
for i0=1:size(idx_val,1)
    formatSpec = msg + "%s\n";
    msg0 =  sprintf("%d ",idx_val(i0,:));
    msg = sprintf(formatSpec,msg0);
end
conf_detail.a4idx_val = msg;

fields = fieldnames(conf_detail);
text_lines = cell(numel(fieldnames(conf_detail)), 1);
for i = 1:numel(fields)
    text_lines{i} = conf_detail.(fields{i});
end

fid = fopen(file_name, 'w');
fprintf(fid, '%s\n', text_lines{:});
fclose(fid);
