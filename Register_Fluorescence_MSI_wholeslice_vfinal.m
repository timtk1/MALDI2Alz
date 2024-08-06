%% Load the MSI data

Animal_3_5xFAD_s1_norm = NegativeDataNorm{5,1};
Animal_3_5xFAD_s2_norm = NegativeDataNorm{6,1};
Animal_1_5xFAD_s1_norm = NegativeDataNorm{1,1};
Animal_1_5xFAD_s2_norm = NegativeDataNorm{1,1};
Animal_2_5xFAD_s2_norm = NegativeDataNorm{4,1};
Animal_4_wt_s1_norm = NegativeDataNorm{7,1};
Animal_4_wt_s2_norm = NegativeDataNorm{8,1};
Animal_5_wt_s2_norm = NegativeDataNorm{10,1};
Animal_6_wt_s1_norm = NegativeDataNorm{11,1};
%Animal_3_5xFAD_s1 = load('Animal_3_5xFAD_s1.mat')
%MSI_data_5xfad = Animal_3_5xFAD_s2_norm;
MSI_data_5xfad = Animal_1_5xFAD_s1_norm;
MSI_data_wt = Animal_4_wt_s1_norm;
load('C:\Projects\AD Effort\Colormaps\Colormaps (5)\Colormaps\viridis')
load('S:\Mar - Imaging - M2 - DHAP\Negative Mode Data\mz_bins_use_neg.mat')

msic = MSI_data_5xfad(:,:,3770); %index m/z channel
msin = msic ./ max(msic(:));


figure(2); imagesc(msin)
colormap(gca,viridis)
xticks([]); yticks([])
clim([0 .8])
colorbar


msic = MSI_data_wt(:,:,3770); %index m/z channel
msin = msic ./ max(msic(:));


figure(3); imagesc(msin)
colormap(gca,viridis)
xticks([]); yticks([])
clim([0 .8])
colorbar

%%
[RGBImage_5xfad,grayscale_image_5xfad] = PCA_MSI(MSI_data_5xfad,[1,2,6]);

%%
[RGBImage_wt,grayscale_image_wt] = PCA_MSI(MSI_data_wt,[1,3,4])


%%

Merged_image = imread("S:\Mar - Imaging - M2 - DHAP\ThioS-Stained Images\May 2\Animal 1 s1 edited\Slide28_Animal1_stitch_merged.png");
Thios_channel = imread("S:\Mar - Imaging - M2 - DHAP\ThioS-Stained Images\May 2\Animal 1 s1 edited\Slide28_Animal1_stitch_channels__Thioflavin S.png");

%change for specific section
% Merged_image = imread("S:\Mar - Imaging - M2 - DHAP\ThioS-Stained Images\May 2\Animal 2 5xfad s2 edited\Slide6_Animal2_merged_.png");
% Thios_channel = imread("S:\Mar - Imaging - M2 - DHAP\ThioS-Stained Images\May 2\Animal 2 5xfad s2 edited\Slide6_Animal2_stitch_channels_Thioflavin S.png");


Merged_imageg = rgb2gray(Merged_image); %convert to grayscale
Thios_channelg = rgb2gray(Thios_channel);

% Merged_image_df = imresize(ThioS_imageg, 0.25); %downsample fluorescence image to increase speed


figure;
imagesc(Merged_imageg)
set(gcf,'position',[213.0000  105.0000  928.0000  651.2000])
xticks([]); yticks([]);
colormap jet
%%

%resive fluorescence images to size of msi image

[height_msi, width_msi] = size(grayscale_image_5xfad);
Merged_imageg_rs = imresize(Merged_imageg, size(grayscale_image_5xfad));
Thios_channel_rs = imresize(Thios_channelg, size(grayscale_image_5xfad));
Thios_channel_rs_ad = imadjust(Thios_channel_rs);

figure; imshow(Merged_imageg_rs)
figure; imshow(Thios_channel_rs_ad)
%%
%Combined_rs = Merged_imageg_rs + Thios_channel_rs_ad*.5;
moving = Merged_imageg_rs ;%./ max(ThioS_imageg_rs(:)) ; %set fluorescence data as moving image
fixed = grayscale_image_5xfad ;%./ max(msi_image(:)); % set MSI as fixed image
% 


figure()
h = imshowpair(moving, fixed)
set(gcf, 'position', [213.0000  105.0000  928.0000  651.2000]);


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
Thios_channel_rs_ad_transform = imwarp(Thios_channel_rs_ad, movingRegisteredAffineWithIC,'OutputView',imref2d(size(fixed))); %register ThioS and MSI data


figure();
imshowpair(ThioS_imageg_rs_transformed,fixed)
set(gcf,'position',[213.0000  105.0000  928.0000  651.2000])


%%

thresholdValue = graythresh(Thios_channel_rs_ad_transform); 

Thios_binary = imbinarize(Thios_channel_rs_ad_transform,thresholdValue);%, 'adaptive',Sensitivity=0.000001); 
%Thios_binary = Thios_channel_rs_ad_transform > 0.5;
figure; imagesc(Thios_binary)
figure()

dilation_element = strel('square', 2); 
Thios_binary_mask_transform_dilated = imdilate(Thios_binary, dilation_element);


%Thios_binary_transform([1:40 123:end],:)=0;
imshowpair(MSI_data_5xfad(:,:,3770),Thios_binary_mask_transform_dilated)
%%

Thiosgbin = Thios_binary_mask_transform_dilated;%imbinarize(Thios_binary_mask_transform_dilated); % binarize i geuss
num_test_pixels = nnz(Thiosgbin);

control_dilation = strel('square', 5);
Thios_dilated_neighbors = imdilate(Thiosgbin, control_dilation);
Thios_control_mask = Thios_dilated_neighbors & ~Thiosgbin;

control_indices = find(Thios_control_mask);
control_indices = control_indices(randperm(length(control_indices), num_test_pixels));

Thios_control_mask_selected = false(size(Thios_control_mask));
Thios_control_mask_selected(control_indices) = true;

figure;
imshowpair(MSI_data_5xfad(:,:,3770),Thios_control_mask_selected)
%%



%%
Thios_binary_transform_rep = repmat(Thios_binary_mask_transform_dilated,[1 1 4031]);
Control_binary_transform_rep = repmat(Thios_control_mask_selected,[1 1 4031]);


MSI_data_plaque_pixels = MSI_data_5xfad.*Thios_binary_transform_rep ;
MSI_data_plaque_pixels2d = reshape(MSI_data_plaque_pixels, [size(MSI_data_plaque_pixels,1)*size(MSI_data_plaque_pixels,2) size(MSI_data_plaque_pixels,3)]   );
meanspec = mean(MSI_data_plaque_pixels2d,1);
Animal_1_5xFAD_s1_plaque_pixels = MSI_data_plaque_pixels;

%%

MSI_data_control_pixels = MSI_data_5xfad.*Control_binary_transform_rep ;
MSI_data_control_pixels2d = reshape(MSI_data_control_pixels, [size(MSI_data_control_pixels,1)*size(MSI_data_control_pixels,2) size(MSI_data_control_pixels,3)]   );

nonZeroRowsplaque = any(MSI_data_plaque_pixels2d ~= 0, 2);
MSI_data_plaque_pixels2d_filtered = MSI_data_plaque_pixels2d(nonZeroRowsplaque, :);

nonZeroRowscontrol = any(MSI_data_control_pixels2d ~= 0, 2);
MSI_data_control_pixels2d_filtered = MSI_data_control_pixels2d(nonZeroRowscontrol,:);

mean_plaque_pixels = mean( MSI_data_plaque_pixels2d_filtered,1);
mean_control_pixels = mean( MSI_data_control_pixels2d_filtered,1);

subspec = mean_plaque_pixels - mean_control_pixels;
figure; stem(mz_bins_use_neg,subspec,'marker','none')

%%
% %%
% 
% % Get the number of original test pixels
% Thiosgbin = imbinarize(Thios_grayscale_dilated); % binarize i geuss
% num_test_pixels = nnz(Thiosgbin);
% 
% control_dilation = strel('square', 5);
% Thios_dilated_neighbors = imdilate(Thiosgbin, control_dilation);
% Thios_control_mask = Thios_dilated_neighbors & ~Thiosgbin;
% 
% control_indices = find(Thios_control_mask);
% control_indices = control_indices(randperm(length(control_indices), num_test_pixels));
% 
% Thios_control_mask_selected = false(size(Thios_control_mask));
% Thios_control_mask_selected(control_indices) = true;
% 
% figure;
% imshowpair(MSI_data(:, :, 3770), Thios_control_mask_selected);
% title('Randomly Selected Control Pixels Adjacent to Test Pixels');
% 
% 
% %%
% MSI_data_wt = NegativeDataNorm{12,1};
% %MSI_data_wt = MSI_data_wt(:,1:217,:);
% MSI_data_plaque_pixelswt = MSI_data_wt.*Thios_binary_transform_rep ;
% MSI_data_plaque_pixels2dwt = reshape(MSI_data_plaque_pixelswt, [size(MSI_data_plaque_pixelswt,1)*size(MSI_data_plaque_pixelswt,2) size(MSI_data_plaque_pixelswt,3)]   );
% meanspecwt = mean(MSI_data_plaque_pixels2dwt,1);
% Animal_6_wt_s2_plaque_pixels = MSI_data_plaque_pixels2dwt;
% %%
% 
% figure;
% plot(mz_bins_use_neg,meanspec);
% hold on;
% plot(mz_bins_use_neg,meanspecwt*-1)
% 
% %%
% 
% figure;
% msi_image()
% imagesc(msi_image)