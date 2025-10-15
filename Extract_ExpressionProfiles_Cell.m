%{
this script loads spatial transcriptomic data from the Cell paper (Chen, 2020) data set
parallel processing is used to speed up this script
%}
dataDir = '.\Cell Paper Data'; %replace with location where data were doanloaded
exprFile = fullfile(dataDir, 'GSE152506_logCPM_counts.txt');
metaFile = fullfile(dataDir, 'spot_metadata.tsv');

% define which slices are APP vs WT 
sliceInfo = struct( ...
    'N02_D1', 'APP', 'N03_D2', 'APP', 'N02_C1', 'APP', 'N03_C2', 'APP', ...
    'N04_D1', 'APP', 'N04_E1', 'APP', 'N05_C2', 'APP', 'N05_D2', 'APP', ...
    'N06_D2', 'APP', 'N07_C1', 'APP', ...
    'B02_D1', 'WT', 'B02_E1', 'WT', 'B03_C2', 'WT', 'B03_D2', 'WT', ...
    'B04_D1', 'WT', 'B04_E1', 'WT', 'B05_D2', 'WT', 'B05_E2', 'WT', ...
    'B06_E1', 'WT', 'B07_C2', 'WT' ...
);

sliceNames = fieldnames(sliceInfo);

% ganglioside genes
genes_synthesis = ["UGCG", "B4GALT5", "B4GALT6", "ST3GAL5", "ST8SIA1", ...
                   "B4GALNT1", "B3GALT4", "ST3GAL2", "ST8SIA5", "GAL3ST1"];
genes_degradation = ["NEU1", "NEU3", "NEU4", "GLB1", "HEXA", "HEXB", ...
                     "ASAH1", "CTSA", "PSAP", "GM2A"];



targetGenes = [genes_synthesis, genes_degradation];

fid = fopen(exprFile, 'r');
headerLine = fgetl(fid); fclose(fid);
headerTokens = strsplit(headerLine, ',');
geneNames = strrep(headerTokens(2:end), '"', '');

[isFound, geneIdx] = ismember(upper(targetGenes), upper(geneNames));
missingGenes = targetGenes(~isFound);
if any(~isFound)
    warning('Missing genes: %s', strjoin(missingGenes, ', '));
end
targetGenes = targetGenes(isFound);
geneCols = geneIdx(isFound);

optsMeta = detectImportOptions(metaFile,'FileType','text','Delimiter','\t');
metaT = readtable(metaFile, optsMeta);

if isempty(gcp('nocreate'))
    parpool;
end

numSlices = numel(sliceNames);
exprCell = cell(numSlices, 1);
groupCell = cell(numSlices, 1);
spotCell  = cell(numSlices, 1);
sliceCell = cell(numSlices, 1);

parfor sIdx = 1:numSlices
    thisSlice = sliceNames{sIdx};
    thisGroup = sliceInfo.(thisSlice);

   fid = fopen(exprFile, 'r');
    fgetl(fid); % skip header

    thisExpr = [];
    thisGroupList = {};
    thisSpots = {};
    thisSlices = {};

    while ~feof(fid)
        line = fgetl(fid);
        tokens = strsplit(line, ',');
        if numel(tokens) < 2, continue; end

        spotName = strrep(tokens{1}, '"', '');
        if startsWith(spotName, [thisSlice '__'])
            values = str2double(tokens(1 + geneCols));
            thisExpr(end+1,:) = values; %#ok<AGROW>
            thisGroupList{end+1,1} = thisGroup;
            thisSpots{end+1,1} = spotName;
            thisSlices{end+1,1} = thisSlice;
        end
    end
    fclose(fid);

    exprCell{sIdx} = thisExpr;
    groupCell{sIdx} = thisGroupList;
    spotCell{sIdx} = thisSpots;
    sliceCell{sIdx} = thisSlices;

end

allExpr  = vertcat(exprCell{:});
allGroup = vertcat(groupCell{:});
allSpot  = vertcat(spotCell{:});
allSlice = vertcat(sliceCell{:});

exprT = array2table(allExpr, 'VariableNames', cellstr(targetGenes));
exprT.Slice = allSlice;
exprT.Group = allGroup;
exprT.Spot  = allSpot;

%%
% %prep stats
fcVals_Cell = zeros(numel(targetGenes),1);
pVals_Cell  = zeros(numel(targetGenes),1);

for i = 1:numel(targetGenes)
    g = targetGenes(i);
    valsAPP = exprT{strcmp(exprT.Group,'APP'), g};
    valsWT  = exprT{strcmp(exprT.Group,'WT'), g};
    [~, p] = ttest2(valsAPP, valsWT);
    pVals_Cell(i) = p;
    fcVals_Cell(i) = mean(valsAPP) - mean(valsWT);
end

