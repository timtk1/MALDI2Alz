%% Load the MSI data

Animal_3_5xFAD_s2 = load('S:\Mar - Imaging - M2 - DHAP\Negative Mode Data\Mat files neg mode/Animal_3_5xFAD_s2.mat')
MSI_data = Animal_3_5xFAD_s2.Animal_3_5xFAD_s2;
load('C:\Projects\AD Effort\Colormaps\Colormaps (5)\Colormaps\viridis')


msic = MSI_data(:,:,1); %index m/z channel
msin = msic ./ max(msic(:));


figure(1); imagesc(msin)
colormap(gca,viridis)
xticks([]); yticks([])
clim([0 .4])
colorbar


%%
%Registration and visualization of ThioS-stained fluorescence data and high
%res (5 um spot size) MSI data


%first, generate PCA representation of data


% Assuming Animal3_S18_HR is your hyperspectral image of size [480 x 410 x 22]
[nRows, nCols, nChannels] = size(MSI_data);
MSI_data2d = reshape(MSI_data, nRows * nCols, nChannels);


% Perform PCA on the 2D data
[coeff, score, latent, tsquared, explained] = pca(MSI_data2d);

% Extract the first principal component scores
firstPC = score(:, 1);

% Reshape the first principal component to the original image size
firstPCImage = reshape(firstPC, nRows, nCols);

% Plot the first principal component
figure;
imagesc(firstPCImage); % Use imagesc to scale the colormap
colormap('gray');      % Set colormap to gray for better visualization
colorbar;              % Show the colorbar
title('First Principal Component');
axis image;            % Set axes to image mode

%%

% Extract the first three principal components
firstPC = score(:, 1);
secondPC = score(:, 2);
thirdPC = score(:, 3);

% Reshape the components to the original image size
firstPCImage = reshape(firstPC, nRows, nCols);
secondPCImage = reshape(secondPC, nRows, nCols);
thirdPCImage = reshape(thirdPC, nRows, nCols);

% Normalize each PC image to [0, 1] for display
firstPCImageNorm = (firstPCImage - min(firstPCImage(:))) / (max(firstPCImage(:)) - min(firstPCImage(:)));
secondPCImageNorm = (secondPCImage - min(secondPCImage(:))) / (max(secondPCImage(:)) - min(secondPCImage(:)));
thirdPCImageNorm = (thirdPCImage - min(thirdPCImage(:))) / (max(thirdPCImage(:)) - min(thirdPCImage(:)));

% Create the RGB image
RGBImage = cat(3, firstPCImageNorm, secondPCImageNorm, thirdPCImageNorm);

% Display the RGB image
figure;
imshow(RGBImage);
xticks([]);  yticks([]);
title('RGB Image top 3 PCs');
axis on;  % Optional: to show axes
set(gcf,'position',[213.0000  105.0000  928.0000  651.2000])
%%
msi_image = rgb2gray(RGBImage);

figure
imshow(grayscale_image)
set(gcf,'position',[213.0000  105.0000  928.0000  651.2000])

%%

%ThioS_image = imread('S:\Mar - Imaging - M2 - DHAP\ThioS-Stained Images\May 2\ThioS_Slide27_Animal3_stitched_c1_bckgrnd.png');
ThioS_image = imread('Animal_3_5xFAD_s2_ThioS.png');

ThioS_imageg = rgb2gray(ThioS_image); %convert to grayscale
ThioS_imageg_ds = imresize(ThioS_imageg, 0.25); %downsample fluorescence image to increase speed



figure;
imshow(ThioS_imageg_ds)
%%

%resive fluorescence images to size of msi image
[height_msi, width_msi] = size(msi_image);
ThioS_imageg_rs = imresize(ThioS_imageg_ds, size(msi_image));

figure; imshow(ThioS_imageg_rs)

%%
moving = ThioS_imageg_rs ./ max(ThioS_imageg_rs(:)) ; %set fluorescence data as moving image
fixed = msi_image./ max(msi_image(:)); % set MSI as fixed image
% 
figure()
imshowpair(moving,fixed)%,"Scaling","independent")
set(gcf,'position',[213.0000  105.0000  928.0000  651.2000])
%% select control points for affine transformation

cpselect(moving,fixed);
%%

% run affine transformation
tformPts = fitgeotrans(movingPoints, fixedPoints, 'affine'); 

[optimizer,metric] = imregconfig('multimodal');
optimizer.MaximumIterations = 500;
optimizer.InitialRadius = 1e-5;
optimizer.Epsilon = 1e-5;

movingRegisteredAffineWithIC = imregtform(moving,fixed,'affine',optimizer,metric,...
    'InitialTransformation',tformPts,'PyramidLevels',3);

ThioS_imageg_rs_transformed = imwarp(moving, movingRegisteredAffineWithIC,'OutputView',imref2d(size(fixed))); %register ThioS and MSI data
figure();
imshowpair(ThioS_imageg_rs_transformed,fixed)
set(gcf,'position',[213.0000  105.0000  928.0000  651.2000])

%%
% %% BLOB DETECTION
% 
% % Assuming ThioS_imageg_ds is your downscaled grayscale image
% % Load and convert the image to grayscale if not already done
% % Assuming ThioS_imageg_ds is your downscaled grayscale image
% 
% % Step 1: Apply the threshold to create a binary image
% thresholdValue = 65;  % Adjust the threshold based on your image's histogram
% binaryImage = moving_transformed > thresholdValue;
% 
% % Step 2: Remove large objects (more than 100 pixels)
% largeBlobMask = bwareaopen(binaryImage, 100);
% finalMask = binaryImage & ~largeBlobMask;
% 
% % Display results
% figure;
% subplot(1, 2, 1);
% imshow(moving_transformed);
% title('Original Image');
% 
% subplot(1, 2, 2);
% imshow(finalMask);
% title('Final Mask');
% 
% %%
% 
% % Assuming finalMask is the mask you created and want to resize
% % Assuming msi_image is already loaded and its size is what you want to match
% 
% % Get the si% Resize the mask using 'nearest' interpolation
% ThioS_imageg_rs_mask = imresize(finalMask, [height_msi, width_msi], 'nearest');
% 
% % Optionally try 'bilinear' or 'bicubic' to see if they might yield better results in this specific case
% % ThioS_imageg_rs_mask = imresize(finalMask, [height_msi, width_msi], 'bilinear');
% 
% % After resizing, apply morphological dilation to recover some lost features
% se = strel('disk', 1);  % Adjust the size based on the extent of detail loss
% ThioS_imageg_rs_mask = imdilate(ThioS_imageg_rs_mask, se);
% 
% % Convert any non-zero values back to logical to ensure it remains a binary mask
% ThioS_imageg_rs_mask = ThioS_imageg_rs_mask > 0;
% 
% % Display the resized and processed mask
% figure;
% imshow(ThioS_imageg_rs_mask);
% title('Resized and Processed Mask');

%%


% Read the fluorescence image
% Read the fluorescence image


% Manually set the threshold value
manualThreshold = 120; % Adjust this value based on your image
binaryMask = ThioS_imageg_ds > manualThreshold;

% Perform blob analysis to identify connected regions
blobAnalysis = regionprops(binaryMask, 'Centroid', 'PixelIdxList','Area');

% Display the original image with detected spots
figure;
imshow(ThioS_imageg_ds);
hold on;
for i = 1:length(blobAnalysis)
    if blobAnalysis(i).Area <= 100
        centroid = blobAnalysis(i).Centroid;
        plot(centroid(1), centroid(2), 'r.'); % Mark the centroid with a red asterisk
    end
end
hold off;
set(gcf,'position',[213.0000  105.0000  928.0000  651.2000])

binaryMask = false(size(binaryMask));  % Initialize a new binary mask of the same size
% Iterate through each blob
for i = 1:length(blobAnalysis)
    % Check the area of each blob
    if blobAnalysis(i).Area <= 100
        % If the blob is smaller than or equal to 100 pixels, include it in the new binary mask
        binaryMask(blobAnalysis(i).PixelIdxList) = true;
    end
end


hold off;

% Display the binary mask
figure;
imshow(binaryMask);

%%

%%

allCoordinates = [];
% Extract pixel coordinates for each detected spot
for i = 1:length(blobAnalysis)
    pixelIndices = blobAnalysis(i).PixelIdxList;
    [rows, cols] = ind2sub(size(binaryMask), pixelIndices);
    %disp(['Spot ', num2str(i), ' - Pixel Coordinates:']);
    %disp([rows, cols]);

    % Extract pixel coordinates for the current spot
    pixelIndices = blobAnalysis(i).PixelIdxList;
    [rows, cols] = ind2sub(size(binaryMask), pixelIndices);
    
    % Append coordinates to the array
    coordinates = [cols, rows]; % [x, y]
    allCoordinates = [allCoordinates; coordinates];
    
end

%%

figure;
tl = tiledlayout(2,2)
nexttile
imagesc(ThioS_imageg_ds);
xticks([]); yticks([])
nexttile
imagesc(binaryMask)
xticks([]); yticks([])
nexttile
imagesc(ThioS_imageg_ds(330:460,1170:1400))
xticks([]); yticks([])
nexttile
imagesc(binaryMask(330:460,1170:1400))
xticks([]); yticks([])
tl.TileSpacing = "compact" 

%%
ThioS_imageg_rs = imresize(ThioS_imageg_ds,size(msi_image));
%%

% Create a copy of the registered image
registeredImage_with_zeros = moving_transformed;

% Convert the 2D coordinates to linear indices
registeredIndices = sub2ind(size(moving_transformed), registeredCoordinates(:, 2), registeredCoordinates(:, 1));

% Set the pixel values at indices to zero
registeredImage_with_zeros(registeredIndices) = 0;
MSI_slice_with_zeros = MSI_slice;
MSI_slice_with_zeros(registeredIndices) = 0;

% Display the registered image with mapped coordinates set to zero
figure;
%imagesc(registeredImage_with_zeros);
xticks([]); yticks([]);
imagesc(MSI_slice_with_zeros)
xticks([]); yticks([]);

%%
% Original image size
[originalHeight, originalWidth] = size(ThioS_imageg_ds);

% Resized image size (assuming msi_image is already loaded)
[newHeight, newWidth] = size(msi_image);

% Calculate scaling factors for the coordinates
scaleX = newWidth / originalWidth;
scaleY = newHeight / originalHeight;

% Resize the original image
ThioS_imageg_rs = imresize(ThioS_imageg_ds, [newHeight, newWidth]);

% Adjust the coordinates of the blobs to the new image size
adjustedCoordinates = allCoordinates * diag([scaleX, scaleY]);

% Extract and map the coordinates to the new size
for i = 1:size(adjustedCoordinates, 1)
    % Ensure the coordinates are within the image bounds
    adjustedCoordinates(i, 1) = min(max(1, round(adjustedCoordinates(i, 1))), newWidth);
    adjustedCoordinates(i, 2) = min(max(1, round(adjustedCoordinates(i, 2))), newHeight);
end

% Create a copy of the resized image to display the adjusted blobs
resizedImage_with_blobs = ThioS_imageg_rs;

% Convert the adjusted 2D coordinates to linear indices in the resized image
adjustedIndices = sub2ind(size(resizedImage_with_blobs), adjustedCoordinates(:, 2), adjustedCoordinates(:, 1));

% Set the pixel values at these indices to a visible value (e.g., maximum intensity to stand out)
resizedImage_with_blobs(adjustedIndices) = max(resizedImage_with_blobs(:));

% Display the resized image with mapped coordinates highlighted
figure;
imagesc(resizedImage_with_blobs);
colormap('gray'); % Change colormap if needed
axis equal;
xticks([]); yticks([]);
title('Resized Image with Mapped Blob Coordinates');

%%% Assuming 'adjustedIndices' contains the indices of blob locations in the resized image

% Select a specific m/z channel for demonstration
mz_channel =3770; % Example, adjust according to your needs

% Extract the specific channel data
specificChannelData = MSI_data(:, :, mz_channel);

% Normalize the channel data for better visualization
normalizedData = mat2gray(specificChannelData);

% Convert the normalized data to RGB (grayscale initially)
rgbImage = repmat(normalizedData, [1, 1, 3]);

% Create a red overlay for blob locations
% Initialize the overlay - using NaNs to indicate no change
redOverlay = nan(size(normalizedData));

% Set the blob locations to red in the overlay
redOverlay(adjustedIndices) = 1; % Full intensity for red channel

% Insert the red overlay into the RGB image
rgbImage(:, :, 1) = max(rgbImage(:, :, 1), redOverlay); % Merge red overlay

% Display the original m/z channel image and the image with enhanced blobs
figure;
subplot(1, 2, 1);
imagesc(specificChannelData);
title(sprintf('Original m/z Channel %d', mz_channel));
colormap('gray'); % Use grayscale for the original image
colorbar;
axis equal;
xticks([]);
yticks([]);

subplot(1, 2, 2);
imshow(rgbImage);
title(sprintf('Blob Highlighted m/z Channel %d', mz_channel));
axis equal;
xticks([]);
yticks([]);

% Annotation to describe the visualization
annotation('textbox', [0.5, 0.01, 0, 0], 'string', 'Blobs highlighted in red', ...
           'HorizontalAlignment', 'center', 'FitBoxToText', 'on', 'EdgeColor', 'none');

%%

% Assuming you already have 'movingRegisteredAffineWithIC' from previous registration
% and 'resizedImage_with_blobs' is ready to be transformed.

% Apply the affine transformation to the blob-highlighted image
resizedImage_with_blobs_transformed = imwarp(resizedImage_with_blobs, movingRegisteredAffineWithIC, 'OutputView', imref2d(size(fixed)));

% Display the registered blob-highlighted image alongside the fixed MSI image for comparison
figure();
imshowpair(resizedImage_with_blobs_transformed, fixed, 'Scaling', 'independent');
title('Registered Blob Image and Fixed MSI Image');
set(gcf, 'position', [213.0000  105.0000  928.0000  651.2000]);

% Optionally, if you want to save or further process the registered image
% imwrite(resizedImage_with_blobs_transformed, 'PathToSaveRegisteredBlobImage.png');


%%

% Ensure the blob mask is binary
blobMask = resizedImage_with_blobs_transformed > 0;

% Dimensions of the MSI data
[numRows, numCols, numChannels] = size(MSI_data);

% Find indices of blob pixels in the mask
blobIndices = find(blobMask); 
numBlobs = length(blobIndices); % Number of blob pixels

% Initialize the matrix to hold the blob data from all channels
blobDataMatrix = zeros(numBlobs, numChannels);

% Extract intensities from each channel at the blob locations
for i = 1:numChannels
    currentChannelData = MSI_data(:, :, i); % Extract data for the current channel
    blobDataMatrix(:, i) = currentChannelData(blobIndices); % Extract blob intensities
end

% Optionally, display the intensities for the first blob across all channels
figure;
plot(blobDataMatrix(1, :)); % Plot the spectral intensity profile of the first blob
xlabel('m/z Channel');
ylabel('Intensity');
title('Intensity Profile of the First Blob Pixel Across All m/z Channels');


%%


Animal_4_wt_s1 = load('S:\Mar - Imaging - M2 - DHAP\Negative Mode Data\Mat files neg mode/Animal_4_wt_s1.mat')
MSI_data_wt = Animal_4_wt_s1.Animal_4_wt_s1;

%%
% Initialize the matrix to hold the blob data from all channels
blobDataMatrix_wt = zeros(numBlobs, numChannels);

% Extract intensities from each channel at the blob locations
for i = 1:numChannels
    currentChannelData = MSI_data_wt(:, :, i); % Extract data for the current channel
    blobDataMatrix_wt(:, i) = currentChannelData(blobIndices); % Extract blob intensities
end


ms5xfad = sum(blobDataMatrix,1);
mswt = sum(blobDataMatrix_wt,1);

figure;
h = stem(mz_bins_use_neg,ms5xfad,'linewidth',1)
hold on
h2 = stem(mz_bins_use_neg,mswt*-1,'linewidth',1)
set(h,'marker','none')
set(h2,'marker','none')

%%
% Number of features
numFeatures = size(blobDataMatrix, 2);

% Initialize arrays to store p-values and fold changes
pValues = zeros(numFeatures, 1);
foldChanges = zeros(numFeatures, 1);

% Assume mz_bins_use_22 is a vector containing m/z values corresponding to each feature
% Ensure mz_bins_use_22 is loaded and accessible here

% Loop through each feature
for i = 1:numFeatures
    % Perform Wilcoxon rank-sum test
    [p, ~, stats] = ranksum(blobDataMatrix(:, i), blobDataMatrix_wt(:, i));
    
    % Store the p-value
    pValues(i) = p;
    
    % Calculate fold change (using median here; replace with mean if more appropriate)
    foldChange = median(blobDataMatrix(:, i)) / median(blobDataMatrix_wt(:, i));
    
    % Store the fold change; log2 transformation
    foldChanges(i) = log2(foldChange);
end
%%
% sizCreate a table with the m/z values, p-values, and fold changes
resultsTable = table(mz_bins_use_neg', pValues, foldChanges, ...
    'VariableNames', {'mz', 'PValue', 'Log2FoldChange'});

% Sort the table by p-value, then by fold change if needed
sortedTable = sortrows(resultsTable, 'PValue');

% Optionally, display the sorted table or write to a file
disp(sortedTable(1:10, :));  % Display the top 10 rows for quick inspection


