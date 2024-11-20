

NegativeDataMarch = {};
NegativeDataMarch{1,1} = Animal_1_5xFAD_s1;
NegativeDataMarch{2,1} = Animal_1_5xFAD_s2;
NegativeDataMarch{3,1} = Animal_2_5xFAD_s1;
NegativeDataMarch{4,1} = Animal_2_5xFAD_s2;
NegativeDataMarch{5,1} = Animal_3_5xFAD_s1;
NegativeDataMarch{6,1} = Animal_3_5xFAD_s2;
NegativeDataMarch{7,1} = Animal_4_wt_s1;
NegativeDataMarch{8,1} = Animal_4_wt_s2;
NegativeDataMarch{9,1} = Animal_5_wt_s1;
NegativeDataMarch{10,1} = Animal_5_wt_s2;
NegativeDataMarch{11,1} = Animal_6_wt_s1;
NegativeDataMarch{12,1} = Animal_6_wt_s2;


%%
% Initialize a variable to hold the combined TIC map
NegativeDataNorm = {};

% Sum the TIC of each data cube to get the combined TIC
for i = 1:12
    data_i = NegativeDataMarch{i,1};
    data_unfold = reshape(data_i,[size(data_i,1).*size(data_i,2)], size(data_i,3)  ) ;
    TIC = squeeze(sum(data_unfold,1));
    lipid_signal = sum(TIC(2860:2920));
    data_norm = data_i ./ lipid_signal;
    NegativeDataNorm{i,1} = data_norm;
end

%%


figure(2)
% Prepare a tiled layout for subplots. Adjust 'flow' to 'grid' if you want a specific layout
tiledlayout(2,6); 

% Iterate through each data cube
for i =1:12
    
    % Extract the current data cube
    currentCube = NegativeDataNorm{i,1};
    
    % Calculate the average spectrum across all X, Y positions for this cube
    % This collapses the X and Y dimensions, leaving the average intensity for each m/z channel
    averageSpectrum = mean(mean(currentCube, 1), 2);
    
    % The result of mean(mean()) is 1x1xZ for m/z channels, so squeeze it to get a 1D array
    averageSpectrum = squeeze(averageSpectrum);
    
    % Plot the average spectrum in a subplot
    nexttile;
    plot(averageSpectrum);
    ylim([0 .1e-4])
   % xlim([3000 4000])
    xlabel('m/z channel index');
    ylabel('Average intensity');


end

% Adjust the layout and add a super title if needed
sgtitle('Average Spectrum of Each Data Cube');




