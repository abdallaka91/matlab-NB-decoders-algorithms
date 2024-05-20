function text_lines  = details_in_lines(ebn0, FER,BER, SER, gen_seq_cnt,K,p, aver_iter, conf_detail)

snr_cnt = length(ebn0);
msg='';
for i0=1:snr_cnt
    formatSpec = msg + "%s";
    msg0 =  sprintf("EbNo = %.3g dB, FER = %d/%d = %g || SER = %d/%d = %g || BER = %d/%d = %g, aver_iter = %g\n" +...
        "------------------------------------------------------------------------------------------------------------------------------------\n",...
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
current_date = datestr(now, 'yyyy_mm_dd');
current_time = datestr(now, 'HH_MM');
file_name = strcat(extractAfter(conf_detail.a1fl_nme, "H matrix : "), '_' ,current_date ,'_' ,current_time, '.txt');

fid = fopen(file_name, 'w');
fprintf(fid, '%s\n', text_lines{:});
fclose(fid)

% clear conf_detail
% conf_detail.a1fl_nme = sprintf("H matrix : %s",H_matrix_mat_fl_nm);
% conf_detail.a2Code = sprintf("N = %d, M = %d, K = %d, GF(%d)",N,M,K,q);
% conf_detail.a3votesW = sprintf('V0 = %f, V1 = %f', v_weights(1) , v_weights(2));
% conf_detail.a4iter = sprintf("H matrix : %d",max_iter);
% conf_detail.a5max_seq = sprintf("max seq generation : %d", max_gen);
% conf_detail.a6fl_nme = sprintf("max error error frame detection: %d",max_err_cnt);
% conf_detail.a7dev_cnt = sprintf("nb of defiation paths : %d", length(dev_lsts_i));
% conf_detail.a8reg_widths = sprintf("regions width (high to low reliable) : {") + sprintf(" %d ", dc11) + sprintf("}");
% formatted_str = '';
% for i = 1:length(di)
%     formatted_str = [formatted_str, sprintf('[%d %d]', di{i}(2:end))];
%     if i < length(di)
%         formatted_str = [formatted_str, ' || '];
%     end
% end
% conf_detail.a9di1_di2 = sprintf("%s\n", formatted_str);