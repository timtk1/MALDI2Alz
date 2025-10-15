
%{
This script pulls in the spatial transcriptomic data from the BMC Genomics
data set extracts expresion profiles for the genes of interest
(genes_synthesis and genes_degradation). To run, download both the
expression data from the paper as well as the region masks (Region Masks
BMC.zip on the IL Databank).
%}

data_dir = '...\Download\ST_dataset'; %change to path of downloaded data
wt_samples = compose("WT.%d", 1:10); % wt animals
tg_samples = compose("TG.%d", 1:10); %5xFAD animals are labeled as TG
all_samples = [wt_samples, tg_samples];

regions = {'Hippocampus', 'CerebralCortex', 'EntorhinalCortex', 'MidBrain'};

% genes to be extracted. these can be modified but make sure case matches
genes_synthesis = ["UGCG", "B4GALT5", "B4GALT6", "ST3GAL5", "ST8SIA1", ...
                   "B4GALNT1", "B3GALT4", "ST3GAL2", "ST8SIA5", "GAL3ST1"];

genes_degradation = ["NEU1", "NEU3", "NEU4", "GLB1", "HEXA", "HEXB", ...
                     "ASAH1", "CTSA", "PSAP","GM2A"];

%concatenate
genes_to_check = [genes_synthesis, genes_degradation];
group_labels = {'WT', 'TG'};

features = readtable(fullfile(data_dir, 'WT.1', 'filtered_feature_bc_matrix', 'features.tsv'), ...
    'FileType', 'text', 'ReadVariableNames', false);
gene_list = string(features.Var2);

%store gene keys in uppercase 
gene_indices = containers.Map();
for g = gene_list'
    gene_indices(upper(g)) = find(strcmpi(gene_list, g), 1);
end

region_expr = struct();
for g = 1:length(genes_to_check)
    gene = genes_to_check(g);
   
    gene_field = matlab.lang.makeValidName(gene);

    for r = 1:length(regions)
        for gr = 1:length(group_labels)
            region_expr.(gene).(regions{r}).(group_labels{gr}) = [];
        end
    end
end

%loops throuh all files , pulls data out. a copule images were orieneted
%differently than the rest and are hardcoded to be corrected
for s = 1:length(all_samples)
    sample_id = all_samples{s};
    base_path = fullfile(data_dir, sample_id);
    sample_key = matlab.lang.makeValidName(sample_id);

    mtx = read_mtx_file(fullfile(base_path, 'filtered_feature_bc_matrix', 'matrix.mtx'));
    T = readtable(fullfile(base_path, 'spatial', 'tissue_positions_list.csv'), 'FileType','text','ReadVariableNames',false);
    barcodes = strtrim(readlines(fullfile(base_path, 'filtered_feature_bc_matrix', 'barcodes.tsv')));
    scalefactors = jsondecode(fileread(fullfile(base_path, 'spatial', 'scalefactors_json.json')));
    scale = scalefactors.tissue_hires_scalef;
    img = imread(fullfile(base_path, 'spatial', 'tissue_hires_image.png'));

    tissue_barcodes = strtrim(string(T.Var1));
    x_raw = T.Var5;
    y_raw = T.Var6;
    [~, idx_expr, idx_spatial] = intersect(barcodes, tissue_barcodes);
    x_scaled = y_raw(idx_spatial) * scale;
    y_scaled = x_raw(idx_spatial) * scale;

    % optional rotation (for viewing images in proper orientation)
    if strcmp(sample_id, 'TG.8')
        [img_h, img_w, ~] = size(img);
        img = flipud(img); % flip vertically
        x_scaled = img_w - x_scaled;
        y_scaled = img_h - y_scaled;

elseif strcmp(sample_id, 'TG.7')
    [img_h, img_w, ~] = size(img);
    x_scaled = img_w - x_scaled;
    x_scaled = x_scaled + 4650; %
    end


    x_pix = round(x_scaled);
    y_pix = round(y_scaled);

    % load masks - these are binary masks that define the brain regions,
    % they can be downloaded from the IL Databank or generate your own
    mask_file = sprintf('region_masks_ST_%s.mat', sample_id);
    if ~isfile(mask_file)
        fprintf(' Mask missing for %s\n', sample_id); continue;
    end
    masks = load(mask_file);

    % in this data set, group IDs begin with WT for wild type animal data and TG
    % for 5xFAD animal data
    if startsWith(sample_id, 'WT')
        group = 'WT';
    elseif startsWith(sample_id, 'TG')
        group = 'TG';
    else
        error("Unrecognized group: %s", sample_id);
    end

    % extract expression per region
    for r = 1:length(regions)
        region = regions{r};
        mask = masks.combinedMasks.(region);

        valid_idx = x_pix > 0 & x_pix <= size(mask,2) & ...
                    y_pix > 0 & y_pix <= size(mask,1);
        x_valid = x_pix(valid_idx);
        y_valid = y_pix(valid_idx);
        in_mask = mask(sub2ind(size(mask), y_valid, x_valid));

        if sum(in_mask) == 0
            continue;
        end

        spot_idx = idx_expr(valid_idx);
        for g = 1:length(genes_to_check)
            gene = genes_to_check(g);
            gene_idx = gene_indices(upper(gene));
            gene_vals = full(mtx(gene_idx, spot_idx))';
            gene_vals_masked = gene_vals(in_mask);
            region_expr.(gene).(region).(group)(end+1) = mean(log1p(gene_vals_masked), 'omitnan');
        end
    end
end

%optional save call
save('region_expr_fixed.mat', 'region_expr', 'genes_synthesis', 'genes_degradation');



%% volcano plot

groups = {'WT', 'TG'};
pVals = zeros(length(genes_to_check),1);
fcVals = zeros(length(genes_to_check),1);

% get animal IDs
wt_samples = compose("WT.%d", 1:10);
tg_samples = compose("TG.%d", 1:10);

for g = 1:length(genes_to_check)
    gene = genes_to_check(g);
    
    % per-animal mean across regions
    mean_expr_WT = zeros(10,1);
    mean_expr_TG = zeros(10,1);
    
    % WT animals
    for i = 1:10
        animal_vals = [];
        for r = 1:length(regions)
            region = regions{r};
            vals = region_expr.(gene).(region).WT;
            animal_vals(end+1) = vals(i);
        end
        mean_expr_WT(i) = mean(animal_vals, 'omitnan');
    end
    
    for i = 1:10
        animal_vals = [];
        for r = 1:length(regions)
            region = regions{r};
            vals = region_expr.(gene).(region).TG;
            animal_vals(end+1) = vals(i);
        end
        mean_expr_TG(i) = mean(animal_vals, 'omitnan');
    end
    
    % run t-test
    [~, p] = ttest2(mean_expr_TG, mean_expr_WT);
    pVals(g) = p;
    
    % calculate log2 fold chance
    fcVals(g) = log2(mean(mean_expr_TG) / mean(mean_expr_WT));
end

figure;
sig_mask = pVals < 0.05;

scatter(fcVals(~sig_mask), -log10(pVals(~sig_mask)), 60, [0.6 0.6 0.6], 'filled'); hold on;
scatter(fcVals(sig_mask), -log10(pVals(sig_mask)), 60, [0.3 0.3 0.9], 'filled');

% top genes are labeled
score = abs(fcVals) .* -log10(pVals); 
[~, sorted_idx] = sort(score, 'descend');
topN = 5;
top_idx = sorted_idx(1:topN);

scatter(fcVals(top_idx), -log10(pVals(top_idx)), 60, 'b', 'filled');
text(fcVals(top_idx), -log10(pVals(top_idx)), cellstr(genes_to_check(top_idx))', ...
    'VerticalAlignment','bottom','HorizontalAlignment','right', ...
    'FontWeight','bold','FontSize',12,'Color','b');

xlabel('Log_2 Fold Change');
ylabel('-log_{10}(p-value)');
title('5xFAD (Lee et al. 2024)','fontweight','normal');

set(gca,'fontsize',14,'linewidth',1,'tickdir','out')
%%


function A = read_mtx_file(filename)
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open the file.');
    end

    header = fgetl(fid);

    line = fgetl(fid);
    while startsWith(line, '%')
        line = fgetl(fid);
    end

    dims = sscanf(line, '%d %d %d');
    rows = dims(1);
    cols = dims(2);
    entries = dims(3);

    data = fscanf(fid, '%d %d %f', [3, entries]);
    fclose(fid);
    
    A = sparse(data(1, :), data(2, :), data(3, :), rows, cols);
end
