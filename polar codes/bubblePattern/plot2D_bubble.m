function plot2D_bubble(matrices, nH, nL, num_levels)

figure; hold on; colormap('hot'); axis off; colorbar;
total_width = sum(1 ./ 2.^(0:num_levels-1));
x_offset = 0;
xticks_pos = [];
xticks_labels = {};
y_ticks_pos = [];
y_ticks_labels = {};

for i = 1:num_levels
    num_matrices = size(matrices{i}, 3);
    scale_factor = 1 / 2^(i-1);
    matrix_wdth = nL * scale_factor;
    pxlwdth = matrix_wdth / nL;
    matrix_hgh = nH * scale_factor;
    pxlhgh = matrix_hgh / nH;

    for j = 1:num_matrices
        x_pos = x_offset;
        y_pos = (j-1) * matrix_hgh;

        img = imagesc([x_pos+pxlwdth/2, x_pos + matrix_wdth-pxlwdth/2],...
            [y_pos+pxlhgh/2, y_pos + matrix_hgh-pxlhgh/2], matrices{i}(:,:,j));

        for k = 0:nL-1
            plot([x_pos + k * pxlwdth, x_pos + k * pxlwdth], ...
                [y_pos, y_pos + matrix_hgh], 'k', 'LineWidth', 0.1);
        end
        for k = 0:nH-1
            plot([x_pos, x_pos + matrix_wdth], ...
                [y_pos + k * pxlhgh, y_pos + k * pxlhgh], 'k', 'LineWidth', 0.1);
        end
    end

    xticks_pos = [xticks_pos, x_offset];
    xticks_labels{end+1} = sprintf('L_{%d}', i);
    x_offset = x_offset + matrix_wdth;
end

xticks(xticks_pos);
xticklabels(xticks_labels);
set(gca, 'TickLabelInterpreter', 'tex');

yticks(y_ticks_pos);
yticklabels(y_ticks_labels);

xlabel('Layer');
ylabel('Cluster');

set(gca, 'YDir', 'reverse'); % Flip figure upside down

hold off;
end
