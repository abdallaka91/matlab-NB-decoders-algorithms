function text_lines  = details_in_lines(ebn0, FER,BER, SER, gen_seq_cnt,K,p, aver_iter, conf_detail, file_name, save_rslt)

snr_cnt = length(ebn0);
msg='\n';
for i0=1:snr_cnt
    formatSpec = msg + "%s";
    msg0 =  sprintf("EbNo = %.3f dB, FER = %4d/%8d = %.3e || SER = %4d/%9d = %.3e || BER = %5d/%11d = %.3e, aver_iter = %.3f\n" +...
        "------------------------------------------------------------------------------------------------------------------------------------------------\n",...
        ebn0(i0), ...
        FER(i0), gen_seq_cnt(i0), FER(i0)/gen_seq_cnt(i0),...
        SER(i0), gen_seq_cnt(i0)*K, SER(i0)/(gen_seq_cnt(i0)*K),...
        BER(i0), gen_seq_cnt(i0)*K*p, BER(i0)/(gen_seq_cnt(i0)*K*p),...
        aver_iter(i0) );
    msg = sprintf(formatSpec,msg0);
end
conf_detail.msg = msg;


fields = fieldnames(conf_detail);
text_lines = cell(numel(fieldnames(conf_detail)), 1);
for i = 1:numel(fields)
    text_lines{i} = conf_detail.(fields{i});
end
if save_rslt
    file_name  =[file_name '.txt'];
    fid = fopen(file_name, 'w');
    fprintf(fid, '%s\n', text_lines{:});
    fclose(fid);
end
