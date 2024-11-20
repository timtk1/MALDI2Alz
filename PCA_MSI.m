function [RGBImage,grayscale_image] = PCA_MSI(MSI_data,components)

%MSI_data is the 3way datacube of MSI data

%components is an array of the 3 principle component the user would like to
%generate the RGB image

%outputs are the RGB image (2D) and grayscale image (2D) for use in image
%registration

%PCA on wild t
[nRows, nCols, nChannels] = size(MSI_data);
MSI_data2d = reshape(MSI_data, nRows * nCols, nChannels);


%PCA
%[coeff, score, latent, tsquared, explained] = pca(MSI_data2d);
[coeff, score] = pca(MSI_data2d);

% pc1 score
firstPC = score(:, 1);

%reshape
firstPCImage = reshape(firstPC, nRows, nCols);

% Plot the first principal component
figure;
imagesc(firstPCImage); % Use imagesc to scale the colormap
colormap('gray');      % Set colormap to gray for better visualization
colorbar;              % Show the colorbar
title('First Principal Component');
axis image;            % Set axes to image mode



% Extract the first three principal components
firstPC = score(:, components(1));
secondPC = score(:, components(2));
thirdPC = score(:, components(3));

% reshape the components to the original image size
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
% tl = tiledlayout(2,1)
% ax = nexttile;
imshow(RGBImage);
xticks([]);  yticks([]);
set(gcf,'position',[213.0000  105.0000  928.0000  651.2000])
% set(ax,'position',[0.3547    0.5875    0.3256    0.3375])
% title('RGB Image top 3 PCs');
% axis on;  % Optional: to show axes
% ax2 = nexttile
% imshow(RGBImage);
% xticks([]);  yticks([]);
% xlim([20 90])
% ylim([36 90])


grayscale_image = rgb2gray(RGBImage);

figure
imshow(grayscale_image)
set(gcf,'position',[213.0000  105.0000  928.0000  651.2000])

end