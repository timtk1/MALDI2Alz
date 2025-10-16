%% Load the MSI data

%available from Illinois Databank - link on Github and in Manuscript
load('mz_bins_use_neg.mat')
load('Animal_1_5xFAD_s1')

MSI_data_5xfad = Animal_1_5xFAD_s1;

msic = MSI_data_5xfad(:,:,3770); %index m/z channel - 3370 is GM3
msin = msic ./ max(msic(:));


figure(2); imagesc(msin)
%colormap(gca,viridis)
xticks([]); yticks([])
clim([0 .8])
colorbar


%%

%{
Run PCA on MSi image. Using the loadings from PCA, we can construct an MSI
image that captures morphology from many MSI 

%}
[RGBImage_5xfad,grayscale_image_5xfad] = PCA_MSI(MSI_data_5xfad,[1,2,6]);

%%
%[RGBImage_wt,grayscale_image_wt] = PCA_MSI(MSI_data_wt,[1,3,4])


%%

% load optical images (download from Illinois databank)
Merged_image = imread('Slide28_Animal1_stitch_merged.png');
Thios_channel = imread('Slide28_Animal1_stitch_channels__Thioflavin S.png');


Merged_imageg = rgb2gray(Merged_image); %convert to grayscale
Thios_channelg = rgb2gray(Thios_channel);


% visualize fluorescence image
figure;
imagesc(Merged_imageg)
set(gcf,'position',[213.0000  105.0000  928.0000  651.2000])
xticks([]); yticks([]);
colormap jet
%%

%resive fluorescence images to size of msi image (fluorescence data is much
%higher spatial resolution)

[height_msi, width_msi] = size(grayscale_image_5xfad);
Merged_imageg_rs = imresize(Merged_imageg, size(grayscale_image_5xfad));
Merged_image_rs_ad = imadjust(Merged_imageg_rs);
Thios_channel_rs = imresize(Thios_channelg, size(grayscale_image_5xfad));
Thios_channel_rs_ad = imadjust(Thios_channel_rs);

figure; imshow(Merged_image_rs_ad)
figure; imshow(Thios_channel_rs_ad)
%%

%setup for image registration

% define fluorescence as moving image, MSI as fixed image. 
moving = Merged_imageg_rs ;%./ max(ThioS_imageg_rs(:)) ; %set fluorescence data as moving image
fixed = grayscale_image_5xfad ;%./ max(msi_image(:)); % set MSI as fixed image
% 

%visualize merged image before registration
figure()
h = imshowpair(moving, fixed)
set(gcf, 'position', [213.0000  105.0000  928.0000  651.2000]);


%% select control points for affine transformation

%{
Select control points. find points in the optical image that match to the
morphological region in the MSI data. we recommend picking points that
include the exterior of the tissue slice, as well as interior morphology
(e.g., hippocampus). After selecting 6-9 control points, run file->export
points to workspace. 
%}

cpselect(moving,fixed);
%%

%{
This block of code runs to affine transformation to match the the
fluorescence image (ThioS channel) to the MSI data
%}

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

%{
This block of code locates all plaques. Image dilation is used to increase
the region by a few pixels to ensure the entire plaque environment is
captured.
%}

thresholdValue = graythresh(Thios_channel_rs_ad_transform); 

Thios_binary = imbinarize(Thios_channel_rs_ad_transform,thresholdValue);%, 'adaptive',Sensitivity=0.000001); 
Thios_binary = imbinarize(Thios_channel_rs_ad_transform,'adaptive',Sensitivity=0.000001); 
figure; imagesc(Thios_binary)
figure()

dilation_element = strel('square', 2); 
Thios_binary_mask_transform_dilated = imdilate(Thios_binary, dilation_element);


%Thios_binary_transform([1:40 123:end],:)=0;
imshowpair(MSI_data_5xfad(:,:,3770),Thios_binary_mask_transform_dilated)
%%

%{
Next, to select control pixels, regions surrounding the plaques are
selected using dilation.
%}

Thiosgbin = Thios_binary_mask_transform_dilated;%imbinarize(Thios_binary_mask_transform_dilated); % binarize
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

%{
Reshaping data
%}
Thios_binary_transform_rep = repmat(Thios_binary_mask_transform_dilated,[1 1 4031]);
Control_binary_transform_rep = repmat(Thios_control_mask_selected,[1 1 4031]);


MSI_data_plaque_pixels = MSI_data_5xfad.*Thios_binary_transform_rep ;
MSI_data_plaque_pixels2d = reshape(MSI_data_plaque_pixels, [size(MSI_data_plaque_pixels,1)*size(MSI_data_plaque_pixels,2) size(MSI_data_plaque_pixels,3)]   );
meanspec = mean(MSI_data_plaque_pixels2d,1);
Animal_1_5xFAD_s1_plaque_pixels = MSI_data_plaque_pixels;

%%

%reshaping control pixels, removing zero-containing rows

MSI_data_control_pixels = MSI_data_5xfad.*Control_binary_transform_rep ;
MSI_data_control_pixels2d = reshape(MSI_data_control_pixels, [size(MSI_data_control_pixels,1)*size(MSI_data_control_pixels,2) size(MSI_data_control_pixels,3)]   );

nonZeroRowsplaque = any(MSI_data_plaque_pixels2d ~= 0, 2);
MSI_data_plaque_pixels2d_filtered = MSI_data_plaque_pixels2d(nonZeroRowsplaque, :);

nonZeroRowscontrol = any(MSI_data_control_pixels2d ~= 0, 2);
MSI_data_control_pixels2d_filtered = MSI_data_control_pixels2d(nonZeroRowscontrol,:);

mean_plaque_pixels = mean( MSI_data_plaque_pixels2d_filtered,1);
mean_control_pixels = mean( MSI_data_control_pixels2d_filtered,1);
