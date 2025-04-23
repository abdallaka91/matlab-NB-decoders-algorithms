clear

rootDir = 'C:\Users\Abdallah Abdallah\Documents\Personals\matlab-NB-decoders-algorithms\polar codes\Genie_Aided_Coefs_Optim\GF256_ccsk_reliab_seq_H2L';
N_list = [64, 32, 16, 8, 4, 2];

for N = N_list
    folderPath = fullfile(rootDir, ['N' num2str(N)]);
    files = dir(fullfile(folderPath, ['mat_N' num2str(N) '_GF256_SNR*.txt']));

    for k = 1:length(files)
        filePath = fullfile(folderPath, files(k).name);

        fid = fopen(filePath, 'r');
        rawLines = {};
        while ~feof(fid)
            line = fgetl(fid);
            rawLines{end+1} = line;
        end
        fclose(fid);


        blocks = {};
        currentBlock = {};
        for i = 1:length(rawLines)
            if isempty(rawLines{i})
                if ~isempty(currentBlock)
                    blocks{end+1} = currentBlock;
                    currentBlock = {};
                end
                blocks{end+1} = {''};
            else
                currentBlock{end+1} = rawLines{i};
            end
        end
        if ~isempty(currentBlock)
            blocks{end+1} = currentBlock;
        end

        changed = 0;
        for i = 1:length(blocks)
            if ~isempty(blocks{i}) && ~all(cellfun(@isempty, blocks{i}))
                numLine = str2num(blocks{i}{1});
                if ~isempty(numLine)
                    blocks{i}{1} = num2str(fliplr(numLine));
                    changed = changed + 1;
                end
                if changed == 4
                    break;
                end
            end
        end


        for i = 1:length(blocks)
            if any(contains(blocks{i}, "File Format"))
                for j = 1:length(blocks{i})
                    if contains(blocks{i}{j}, '1-')
                        blocks{i}{j} = "1- most to least reliable channels (taken entropies)";
                    elseif contains(blocks{i}{j}, '2-')
                        blocks{i}{j} = "2- most to least reliable channels (taken error probabilities)";
                    elseif contains(blocks{i}{j}, '3-')
                        blocks{i}{j} = "3- channels entropies (flipped)";
                    elseif contains(blocks{i}{j}, '4-')
                        blocks{i}{j} = "4- channels error probabilities (flipped)";
                    end
                end
                break;
            end
        end

        ln1 = 0;
        fid = fopen(filePath, 'w');
        for i = 1:length(blocks)
            for j = 1:length(blocks{i})
                if all(isstrprop(blocks{i}{j}, 'digit') | isspace(blocks{i}{j}) | blocks{i}{j} == '.' | blocks{i}{j} == 'e' | blocks{i}{j} == '-' | blocks{i}{j} == '+')
                    nums = str2num(blocks{i}{j});
                    
                    if ~isempty(nums)
                        ln1=ln1+1;
                        if ln1>=3
                            fprintf(fid, ['%.12f '], nums);
                        else
                            fprintf(fid, ['%d '], nums);
                        end
                        fprintf(fid, '\n');
                        continue;
                    end
                end
                fprintf(fid, '%s\n', blocks{i}{j});
            end
        end
        fclose(fid);
    end
end
