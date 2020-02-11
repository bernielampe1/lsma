function mse = plot_fcls(results, saveas_prefix)

    [~, ~, p] = size(results);

    alpha = char(97:122);
    mlabels = {};
    for i = 1:p
        mlabels{end+1} = sprintf('(%c)', alpha(i));
    end

    % Each row, one endmember will be the main endmember,
    % 1:Alunite, 2:Buddingtonite, 3:Calcite, 4:Kaolinite, 5:Muscovite
    panel_dim = [4, 2, 2, 1, 1];
    num_per_row = sum(panel_dim.^2);
    offset = 30;
    
    % this graph will plot the 26 panel pixels for each material
    % material A will plot 26 panels from row 1
    % material B will plot 26 panels from row 2
    % ...
    
    % loop over the materials (material number is same as row number)
    A = zeros(num_per_row, p);
    for r = 1:5
        k = 1;
        for c = 1:5
            for i = 1:panel_dim(c)
                for j = 1:panel_dim(c)
                    A(k, r) = results(r * offset + i, c * offset + j, r);
                    k = k + 1;
                end
            end
        end
    end

    A = flipud(A);
    figure; hold on;
    for m = 1:p
        stem3(ones(1, num_per_row)*m, 1:num_per_row, A(:,m))
    end
    grid on;
    view(-36, 40);
    title('All Panel Pixels, All Materials');
    xticklabels(mlabels);
    saveas(gcf, sprintf('%s_allmaterials.png', saveas_prefix));
    
    % accuracy
    A(isnan(A)) = 0;
    T = [0.25 * ones(1, p); 0.5 * ones(5, p); ones(20, p)];
    mse = sum((T(:)-A(:)).^2)/(26*5);

    % these graphs will plot the 26 panel pixels for each material
    % material A will plot 130 panels from all rows
    % material B will plot 130 panels from all rows
    % ...
    
    % loop over the materials (material number is same as row number)
    for z = 1:p
        A = zeros(num_per_row, 5);
        for r = 1:5
            k = 1;
            for c = 1:5
                for i = 1:panel_dim(c)
                    for j = 1:panel_dim(c)
                        A(k, r) = results(r * offset + i, c * offset + j, z);
                        k = k + 1;
                    end
                end
            end
        end

        A = flipud(A);
        figure; hold on;
        for m = 1:5
            stem3(ones(1, num_per_row)*m, 1:num_per_row, A(:,m))
        end
        grid on;
        view(-36, 40);
        title(sprintf('Material %s', mlabels{z}));
        xticklabels(mlabels);
        saveas(gcf, sprintf('%s_material_stem_%d.png', saveas_prefix, z));
    end
    
    % plot the image strengths
%    for i = 1:p
%        figure;
%        imagesc(results(:,:,i))
%        colormap('jet')
%        colorbar;
%        title(sprintf('Material %s', mlabels{i}));
%        saveas(gcf, sprintf('%s_material_image_%d.png', saveas_prefix, i));
%    end
end
