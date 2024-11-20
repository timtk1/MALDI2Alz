%load array of m/z bins
load('mz_bins_use_neg.mat')
load('Animal_1_5xFAD_s1.mat')

%load MSI data
%MSI_data = NegativeDataNorm{1,1};
MSI_data = Animal_1_5xFAD_s1;

%select m/z to plot. this only affects the visualization used for drawing
%regional boundaries, any m/z with morphological detail can be used
mzuse = find(abs(mz_bins_use_neg - 885.54) == min(abs(mz_bins_use_neg - 885.54)));
msi_slice = MSI_data(:,:,mzuse);

figure;
imagesc(msi_slice); colormap(gray);

%%

%define brain regions

%for hippocampus, cortex, and EC, draw two ROIs for each region on each
%side of the brain. for midbrain, only 1 region is drawn

regions = {'Hippocampus', 'CerebralCortex', 'EntorhinalCortex', 'MidBrain'};
combinedMasks = struct();

for i = 1:length(regions)
    region = regions{i};

    fprintf('Draw polygons for %s\n', region);  % Inform the user which region to draw

    % Initialize figure
    figure;
    imagesc(msi_slice); colormap(gray); axis image;

    % Draw the first polygon ROI for the current region
    title(sprintf('Draw the first polygon ROI for %s', region));
    h1 = drawpolygon;
    wait(h1); % Wait for the user to finish drawing

    % Create mask for the first ROI
    mask1 = poly2mask(h1.Position(:,1), h1.Position(:,2), size(msi_slice,1), size(msi_slice,2));

    if ~strcmp(region, 'MidBrain')
        % If the region is not MidBrain, draw a second polygon ROI
        title(sprintf('Draw the second polygon ROI for %s', region));
        h2 = drawpolygon;
        wait(h2); % Wait for the user to finish drawing

        % mask for the second ROI
        mask2 = poly2mask(h2.Position(:,1), h2.Position(:,2), size(msi_slice,1), size(msi_slice,2));

        % combine the two masks for this region
        combinedMasks.(region) = mask1 | mask2;
    else
        % only first mask for midbrain
        combinedMasks.(region) = mask1;
    end

    % Optional: Show the combined mask overlay on the original slice for verification
    figure; 
    imshowpair(msi_slice, combinedMasks.(region), 'montage');
    title(sprintf('%s - Combined Mask', region));
end

%%

this section extracts the normalized intensity of each m/z from each brain
region open the struct "normalizedIntensityForRegions". each brain region
consists of 1xn double where n is the array of m/z channels
(targetMzValues)

%specify m/z channels for analysis
targetMzValues = [1179.73,1207.77,1382.81,1410.84,1544.87,1572.8,642.24,670.27,647.43,673.48,480.30,571.24,599.31,619.13,...
               857.41,885.54,862.45,878.53, 888.54,890.63,904.61,906.63];

%find closest values
closestMzIndices = arrayfun(@(x) findClosestMz(x, mz_bins_use_neg), targetMzValues);

normalizedIntensityForRegions = struct();

%iterate
regionNames = fieldnames(combinedMasks)';
for regionName = regionNames
    % mask for current region
    combinedMask = combinedMasks.(regionName{1});
    
    % calculate num pixels
    numPixels = sum(combinedMask(:));
    
    normalizedIntensity = zeros(1, length(targetMzValues));
    
    %iterate through the closest m/z indices
    for i = 1:length(closestMzIndices)
        index = closestMzIndices(i);
        
       
        specificSlice = MSI_data(:, :, index);
        
        %apply mask
        maskedSlice = specificSlice .* combinedMask;
        
        %calculate the total signal intensity for the masked area
        totalIntensity = sum(maskedSlice(:));
        
        % normalize intentisty by pixel number
        normalizedIntensity(i) = totalIntensity / numPixels;
    end
    
    % store the normalized intensity array in the struct for the current region
    normalizedIntensityForRegions.(regionName{1}) = normalizedIntensity;
end


function closestIndex = findClosestMz(targetValue, mzArray)
    [~, closestIndex] = min(abs(mzArray - targetValue));
end
