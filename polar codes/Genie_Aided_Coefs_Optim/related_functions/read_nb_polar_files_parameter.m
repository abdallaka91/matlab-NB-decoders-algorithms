function [p, data] = read_nb_polar_files_parameter(filename)
fid = fopen(filename, 'r');
if fid == -1
    error('File could not be opened');
end

p = 6;
data = struct();
i = 0;

line = fgets(fid);

while ischar(line)
    if contains(line, 'rate_match')
        cleaned_line = strtrim(line);
        tokens = strsplit(cleaned_line, {' ', ':', '[', ']'});
        if length(tokens) > 3
            p = sscanf(tokens{4}, '%d');
        end
    end

    if startsWith(line, 'src_indx_snr_')
        i=i+1;
        SNR = sscanf(line, 'src_indx_snr_%f');
        data(i).SNR = SNR;
        data(i).Pe = [];
        data(i).channel_sorting = [];
        line = fgets(fid);
        channel_sorting = sscanf(line, '%d');
        data(i).channel_sorting = channel_sorting;
        line = fgets(fid);
        Pe = sscanf(line, '%f');
        data(i).Pe = Pe;
    end

    line = fgets(fid);
    if  line ~= -1
        if startsWith(line, 'gf_coef_snr')
            gf_coef = [];
            line = fgets(fid);
            c1=true;
            while ischar(line) && c1
                row_data = sscanf(line, '%d');
                if ~isempty(row_data)
                    if isempty(gf_coef)
                        j=1;
                        N2 = length(row_data);  % First row to determine the number of columns
                        n = log2(2 * N2);
                        gf_coef = zeros(N2, n);  % Preallocate matrix
                    end
                    gf_coef(:,j) = row_data;
                    line = fgets(fid);
                    j=j+1;
                    if j>n
                        c1=false;
                    end
                end
            end
            data(i).gf_coef = gf_coef;
        end
    end
end

fclose(fid);
end
