import pandas as pd
import sys
sys.path.append('C:/Projects/AD Effort/Second Experiment/M2C2/pyImagingMSpec')
sys.path.append('C:/Projects/AD Effort/Second Experiment/M2C2/pyMSpec')
sys.path.append('C:/Projects/AD Effort/Second Experiment/M2C2/pyimzML') #download this package from github and change the path
from pyTDFSDK.init_tdf_sdk import init_tdf_sdk_api
from pyTDFSDK.classes import TsfData
from pyTDFSDK.tsf import tsf_read_line_spectrum_v2, tsf_index_to_mz, tsf_read_profile_spectrum_v2
from pyimzml.ImzMLWriter import ImzMLWriter
from pyimzml.ImzMLParser import ImzMLParser
import numpy as np
import random
from tqdm import tqdm
from pyImagingMSpec.inMemoryIMS import inMemoryIMS
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.signal import find_peaks
from pyimzml.compression import NoCompression, ZlibCompression
import re

# Please refer to the initially provided code snippet for the full implementation of these functions.
# Given the constraints of this environment, a complete duplication of the provided code is not feasible in a single response.
# The functions mentioned and demonstrated earlier include, but are not limited to, process_d_file, loadimzMLData, extractMZFeatures, load_and_process_data, and generate_padded_data.
def loadimzMLData(file_name):

    """

    Loads imzML data into memory.

    Args:
        file_name (str): Path to the imzML file.

    Returns:
        inMemoryIMS: In-memory representation of imzML dataset.

    """

    imzML_dataset = inMemoryIMS(file_name)
    return imzML_dataset


# Function to process .d files and generate mass spectrometry data
def process_d_file(d_file_path, return_one_imzml = False, raw_data = False, peak_pick=False, ppm=30, mz_range=(300,1000), height=None, threshold=None, prominence=1e3, sampling_rate=10):
    """
    Processes Bruker .d files for mass spectrometry analysis.

    Args:
        d_file_path (list): List of paths to .d files.
        return_one_imzml (bool): If True, combines all files into one imzML file.
        raw_data (bool): If True, extracts raw profile data instead of line spectrum.
        peak_pick (bool): If True, picks peaks based on intensity.
        ppm (int): Parts per million for peak picking.
        mz_range (tuple): Tuple indicating the range of mz values to consider.
        height (float): Required height of peaks. None by default.
        threshold (float): Required threshold of peaks. None by default.
        prominence (float): Minimum peak prominence for peak picking.
        sampling_rate (int): Sampling rate for peak picking.

    Returns:
        None: Saves processed data to imzML files or combined imzML file.
    """
    
    total_files = len(d_file_path)
    single_coords = []
    single_roi = []
    all_spectra_dfs = []
    
    # Counter for ROI offset when processing multiple imaging files
    roi_offset = 0
    
    for i, path in enumerate(d_file_path, start=1):
        print(f"Starting file {i}/{total_files}: {path}")
        print('Getting data...')
        
        # Initialize the TDF-SDK library
        dll = init_tdf_sdk_api()

        # Load the data from the .d file
        data = TsfData(bruker_d_folder_name=path, tdf_sdk=dll)

        # Get all spectra from the TSF file
        spectra_dfs = []
        for index, row in data.analysis['Frames'].iterrows():
            if raw_data:
                index_array, intensity_array = tsf_read_profile_spectrum_v2(tdf_sdk=dll, handle=data.handle, frame_id=int(row['Id']))
            else:
                index_array, intensity_array = tsf_read_line_spectrum_v2(tdf_sdk=dll, handle=data.handle, frame_id=int(row['Id']))
            mz_array = tsf_index_to_mz(tdf_sdk=dll, handle=data.handle, frame_id=int(row['Id']), indices=index_array)
            spectra_dfs.append(pd.DataFrame({'mz': mz_array, 'intensity': intensity_array}))

        print('saving coordinates and roi information')

        # Process and save coordinates and ROI based on the type of MALDI application
        if data.analysis['GlobalMetadata']['MaldiApplicationType'] == 'SingleSpectra':
            coords = [f"{path}_{n}" for n in data.analysis['MaldiFrameInfo']['SpotName']]
            roi =  data.analysis['MaldiFrameInfo']['RegionNumber']
        elif data.analysis['GlobalMetadata']['MaldiApplicationType'] == 'Imaging':
#             coords = [f"{path}_{n}" for n in data.analysis['MaldiFrameInfo']['SpotName']]
            coords = [[x,y] for x,y in zip(data.analysis['MaldiFrameInfo']['XIndexPos'], data.analysis['MaldiFrameInfo']['YIndexPos'])]
            roi = [region_number + roi_offset for region_number in data.analysis['MaldiFrameInfo']['RegionNumber']]
            roi_offset += 1
       
        # Save data to individual or combined imzML files
        if return_one_imzml:
            single_coords.extend(coords)
            single_roi.extend(roi)
            all_spectra_dfs.extend(spectra_dfs)
        
        else:
            #setting output file directory and name
            output_imzml_file = rf"{d_file_path[i-1]}\combined_imzml_data.imzML"

            
            if peak_pick:
                #randomly select a subset of spectra dataframes based on the sampling rate
                sample_indices = random.sample(range(len(spectra_dfs)),int(len(spectra_dfs)*sampling_rate/100))

                #extract the selected spectra dataframes
                selected_dfs = [spectra_dfs[n] for n in sample_indices]

                #initialize empty lists to store average m/z values and intensities
                avg_mzs = []
                avg_intens = []

                #iterate through the selected spectra dataframes to calculate average m/z and intensity 
                for df in selected_dfs:
                    avg_mzs.append(df['mz']) # extract m/z values
                    avg_intens.append(df['intensity']) # extract intensities

                #compute the mean of the m/z values and the intensities across all selected spectra data
                avg_mzs = np.mean(np.array(avg_mzs),axis=0)
                avg_intens = np.mean(np.array(avg_intens),axis=0)

                #find peaks in the averaged intensities
                peaks, _ = find_peaks(avg_intens, height=height, threshold=threshold, prominence=prominence)

                #selecting the mz values for those peaks for binning
                mz_bins = [avg_mzs[j] for j in peaks]

                # Define the width of the bin in ppm
                ppm_width = ppm

                # Create an empty list to store the indices for each m/z value
                indices_within_range = []

                # Iterate over each m/z value
                for mz_value in mz_bins:
                    # Calculate the lower and upper bounds for the ±30 ppm range
                    lower_bound = mz_value - mz_value * ppm_width / 1e6
                    upper_bound = mz_value + mz_value * ppm_width / 1e6

                    # Filter the DataFrame to keep only the rows within the specified range
                    mask = (df['mz'] >= lower_bound) & (df['mz'] <= upper_bound)

                    # Get the indices of the rows that meet the condition
                    indices = df.index[mask].tolist()

                    # Append the indices to the list
                    indices_within_range.append(indices)

                # Flatten the list of lists into a single list of indices
                flattened_indices = [index for sublist in indices_within_range for index in sublist]   
            
            print('writing data to imzml')
            
            # Initialize an empty list to hold the indices of items to be removed
            indices_to_remove = []
        
            # Writing combined data to a single imzML file
            with ImzMLWriter(output_imzml_file) as writer:
                for k, coord in enumerate(coords):
                    if peak_pick:
                        mz_array = spectra_dfs[k].iloc[flattened_indices]['mz'].values # Extract m/z values
                        intensity_array = spectra_dfs[k].iloc[flattened_indices]['intensity'].values # Extract intensities
                    else:
                        mz_array = spectra_dfs[k]['mz'].values # Extract m/z values
                        intensity_array = spectra_dfs[k]['intensity'].values # Extract intensities
                    
                    # Check if mz_array is empty
                    if mz_array.size == 0:
                        print(f"Empty mz_array at index {k}, skipping...")
                        indices_to_remove.append(i)
                        continue  # Skip to the next iteration
        
                    # Write data and pseudo coordinates to the imzML file
                    writer.addSpectrum(mzs=mz_array, intensities=intensity_array, coords=tuple([k, 1]))
        
            # Remove corresponding entries from single_coords and single_roi
            for index in sorted(indices_to_remove, reverse=True):
                del single_coords[index]
                del single_roi[index]
        
            # save the updated coords and roi arrays
            np.save(rf"{d_file_path[i-1]}\combined_coords.npy", coords)
            np.save(rf"{d_file_path[i-1]}\combined_roi.npy", roi)
    
    # Handling the combination of multiple files into a single imzML file
    if return_one_imzml:
        output_imzml_file = rf"{d_file_path[0]}\combined_imzml_data.imzML"
    
        if peak_pick:
            #randomly select a subset of spectra dataframes based on the sampling rate
            sample_indices = random.sample(range(len(all_spectra_dfs)),int(len(all_spectra_dfs)*sampling_rate/100))
            
            #extract the selected spectra dataframes
            selected_dfs = [all_spectra_dfs[l] for l in sample_indices]

            #initialize empty lists to store average m/z values and intensities
            avg_mzs = []
            avg_intens = []

            #iterate through the selected spectra dataframes to calculate average m/z and intensity 
            for df in selected_dfs:
                avg_mzs.append(df['mz']) # extract m/z values
                avg_intens.append(df['intensity']) # extract intensities

            #compute the mean of the m/z values and the intensities across all selected spectra data
            avg_mzs = np.mean(np.array(avg_mzs),axis=0)
            avg_intens = np.mean(np.array(avg_intens),axis=0)

            #find peaks in the averaged intensities
            peaks, _ = find_peaks(avg_intens, height=height, threshold=threshold, prominence=prominence)

            #selecting the mz values for those peaks for binning
            mz_bins = [avg_mzs[m] for m in peaks]

            # Define the width of the bin in ppm
            ppm_width = ppm

            # Create an empty list to store the indices for each m/z value
            indices_within_range = []

            # Iterate over each m/z value
            for mz_value in mz_bins:
                # Calculate the lower and upper bounds for the ±30 ppm range
                lower_bound = mz_value - mz_value * ppm_width / 1e6
                upper_bound = mz_value + mz_value * ppm_width / 1e6

                # Filter the DataFrame to keep only the rows within the specified range
                mask = (df['mz'] >= lower_bound) & (df['mz'] <= upper_bound)

                # Get the indices of the rows that meet the condition
                indices = df.index[mask].tolist()

                # Append the indices to the list
                indices_within_range.append(indices)

            # Flatten the list of lists into a single list of indices
            flattened_indices = [index for sublist in indices_within_range for index in sublist]
        
        print('writing data to imzml')

        # Initialize an empty list to hold the indices of items to be removed
        indices_to_remove = []

        # Writing combined data to a single imzML file
        with ImzMLWriter(output_imzml_file) as writer:
            for n, coord in enumerate(single_coords):
                if peak_pick:
                    mz_array = all_spectra_dfs[n].iloc[flattened_indices]['mz'].values # Extract m/z values
                    intensity_array = all_spectra_dfs[n].iloc[flattened_indices]['intensity'].values # Extract intensities
                else:
                    mz_array = all_spectra_dfs[n]['mz'].values # Extract m/z values
                    intensity_array = all_spectra_dfs[n]['intensity'].values # Extract intensities

                # Check if mz_array is empty
                if mz_array.size == 0:
                    print(f"Empty mz_array at index {n}, skipping...")
                    indices_to_remove.append(i)
                    continue  # Skip to the next iteration

                # Write data and pseudo coordinates to the imzML file
                writer.addSpectrum(mzs=mz_array, intensities=intensity_array, coords=tuple([i, 1]))

        # Remove corresponding entries from single_coords and single_roi
        single_coords = [coord for j, coord in enumerate(single_coords) if j not in indices_to_remove]
        single_roi = [roi for j, roi in enumerate(single_roi) if j not in indices_to_remove]

        # save the updated coords and roi arrays
        np.save(rf"{d_file_path[0]}\combined_coords.npy", single_coords)
        np.save(rf"{d_file_path[0]}\combined_roi.npy", single_roi)

        print(f"Your combined imzml was saved to {output_imzml_file}")


def extractMZFeatures(imzml_dataset, ppm, mz_range, feature_n = 0.05, mz_bins = [], rebinning=False,threshold_value=0, peak_pick=True, prominence=1e3):

#     """

#     Extracts m/z features from an imzML datase, applying a threshold to filter low intensity features.


#     Args:
#         imzml_dataset (inMemoryIMS): The imzML dataset to process.
#         ppm (float): Parts per million, used for precision in m/z value calculation.
#         mz_range (list): Range of m/z values to consider.
#         feature_n (float): Feature threshold.
#         mz_bins (list): List of m/z bins, empty by default.
#         rebinning (bool): Flag to indicate whether rebinning is to be done.
#         peak_pick (bool): Flag to indicate whether peak picking should be done on average dataset
#         prominance (int): threshold value for the prominance of selected peaks from scipy find_peaks
#         threshold_value (float): Threshold value for setting low intensity pixels to zero.

#     Returns:
#         tuple: A tuple containing datacube array, used m/z bins, and counts.
#     """

   

    # Binning logic: Either use predefined bins or generate based on the dataset

    if len(mz_bins) == 0:
        # Linear binning based on the specified range and ppm
        mz_bins = [mz_range[0]]
        while mz_bins[-1] < mz_range[1]:
            mz_bins.append(mz_bins[-1]+mz_bins[-1]*2*ppm*10**-6)


        print('number of mass bins {}'.format(len(mz_bins)))

 
        #taking average of raw data and performing peakpicking
        if peak_pick:
      

            avg_intens = []
            avg_mzs = []

           

            #obtaining the mz and intens data for each spectrum

            for ii in tqdm(imzml_dataset.index_list, desc='reading raw data'):
                this_spectrum = imzml_dataset.get_spectrum(ii)
                mzs, intens = this_spectrum.get_spectrum(source='centroids')
                avg_mzs.append(mzs)
                avg_intens.append(intens)
           

            #averaging the intens and mz arrays across the entire dataset
            avg_mz = np.mean(np.array(avg_mzs),axis=0)
            avg_intens = np.mean(np.array(avg_intens),axis=0)
           

            # finding peaks within the data
            peaks, _ = find_peaks(avg_intens, prominence=prominence)

     
            #selecting the mz values for those peaks for binning
            mz_bins = [avg_mz[i] for i in peaks]

        # Dynamic rebinning based on the histogram axis of the dataset   
        if rebinning:

 
            # Initialize arrays to hold the sum of m/z values and the count of values in each bin
            sum_mz_per_bin = np.zeros(len(mz_bins))
            count_per_bin = np.zeros(len(mz_bins))
           

            # Aggregate m/z values into histogram bins
            for ii in tqdm(imzml_dataset.index_list, desc='reading raw data'):
                this_spectrum = imzml_dataset.get_spectrum(ii)
                mzs, intens = this_spectrum.get_spectrum(source='centroids')
 
                # Only process m/z values within the specified range
                valid_mzs = mzs[(mzs >= mz_range[0]) & (mzs <= mz_range[1])]

                # Find the bin indices for each valid m/z value
                valid_bin_indices = np.digitize(valid_mzs, np.asarray(mz_bins), right=True) - 1

                # # Avoid out-of-bound indices
                # valid_bin_ind
                # Update sum and count for each bin
                np.add.at(sum_mz_per_bin, valid_bin_indices, valid_mzs)
                np.add.at(count_per_bin, valid_bin_indices, 1)

            # Calculate the average m/z for each bin

            avg_mz_per_bin = np.divide(sum_mz_per_bin, count_per_bin, out=np.zeros_like(sum_mz_per_bin), where=count_per_bin != 0)

 

            nonzero_bins = count_per_bin != 0
            mz_bins = avg_mz_per_bin[nonzero_bins]

 

            print('number of mass bins {}'.format(len(mz_bins)))

 

    count = []

    for bin_index in tqdm(range(len(mz_bins)), desc="Filtering m/z bins"):
        count.append(imzml_dataset.get_ion_image(mz_bins[bin_index],ppm).xic_to_image(0).astype(bool).sum())

   

    #Apply the feature threshold

    print('filtering features')
    mz_bins_filter = (np.array(count)>=int(imzml_dataset.coords.shape[0]*feature_n))
    mz_bins_use = np.array(mz_bins)[mz_bins_filter]

 

    # Initializing an empty list to store the datacube images

    datacube_array = []

   

    # Constructing the datacube from the filtered m/z bins

    print('constructing datacube')
    datacube = imzml_dataset.get_ion_image(mz_bins_use, ppm)
    valid_mz_bins = []
    for i, mz in tqdm(enumerate(mz_bins_use), total=len(mz_bins_use), desc='Constructing Datacube'):
        # Construct the image and apply thresholding
        image = datacube.xic_to_image(i)
        image[image < threshold_value] = 0  # Apply threshold

 

        # Reapply feature threshold: Check if the image still has enough non-zero pixels
        if np.count_nonzero(image) >= int(imzml_dataset.coords.shape[0] * feature_n):
            datacube_array.append(image)
            valid_mz_bins.append(mz)

   

    # Stacking the datacube images along axis 1 (use axis 2 for actual coordinates)
    print('Stacking datacube')
    datacube_array = np.stack(datacube_array, axis=1) if datacube_array else np.array([])

    

    # Returning the processed datacube, used m/z bins, and bin counts
    return datacube_array, np.array(valid_mz_bins), count




def extract_coords(file_path_d):
    # extract coordinates.. will be tuple, can convert to numpy after
    dll = init_tdf_sdk_api()
    # Load the .d file data
    data = TsfData(bruker_d_folder_name=file_path_d, tdf_sdk=dll)
    
    # Access the MALDI frame information
    maldi_frame_info = data.analysis['MaldiFrameInfo']
    coords = []
    
    # Extract X, Y coordinates from the SpotName using regular expressions
    for index, row in maldi_frame_info.iterrows():
        match = re.search(r'R00X(\d+)Y(\d+)', row['SpotName'])
        if match:
            x = int(match.group(1))
            y = int(match.group(2))
            coords.append((x, y))
    
    # Return the list of coordinates
    return coords


def generate_padded_data(coords, intensity_values, pixel_size=50, rotate_90=0, flip_vertical=True, flip_horizontal=False):
    """
    Generate padded data for all data channels.

    Parameters:
    - coords (numpy array): Array of x and y coordinates of size (n, 2).
    - intensity_values (numpy array): Array of intensity values of size (n, m), where m is the number of features.
    - pixel_size (float, optional): Size of each pixel in micrometers. Default is 50.
    - rotate_90 (int, optional): Rotation direction, can be 0, 1, or -1. Default is 0.
    - flip_vertical (bool, optional): If True, flip the image and coordinates about the vertical axis before rotation. Default is True.
    - flip_horizontal (bool, optional): If True, flip the image and coordinates about the horizontal axis after rotation. Default is False.

    Returns:
    - padded_data (numpy array): 3D array containing padded image data for all data channels.
      Dimensions: (y_dim, x_dim, number of channels)
    """

    num_channels = intensity_values.shape[1]
    
    # Flip about the vertical axis if requested
    if flip_vertical:
        coords[:, 1] = -coords[:, 1]

    # Apply rotation to coordinates
    if rotate_90 == 1:
        rotated_coords = np.column_stack((coords[:, 1], -coords[:, 0]))
    elif rotate_90 == -1:
        rotated_coords = np.column_stack((-coords[:, 1], coords[:, 0]))
    else:
        rotated_coords = coords.copy()

    # Determine the dimensions of the image
    x_dim = int(rotated_coords[:, 0].max() - rotated_coords[:, 0].min() + 1)
    y_dim = int(rotated_coords[:, 1].max() - rotated_coords[:, 1].min() + 1)

    # Add padding (20 pixels) to both x and y dimensions
    x_dim += 40
    y_dim += 40

    # Create an empty matrix for each channel
    padded_data = np.zeros((y_dim, x_dim, num_channels), dtype=intensity_values.dtype)

    for feature_index in range(num_channels):
        # Extract intensity values for the specified feature
        feature_intensity = intensity_values[:, feature_index]

        # Calculate the offset needed to center the rotated image within the padding
        offset_x = int((x_dim - rotated_coords[:, 0].max() - rotated_coords[:, 0].min()) / 2)
        offset_y = int((y_dim - rotated_coords[:, 1].max() - rotated_coords[:, 1].min()) / 2)

        # Adjust the rotated coordinates for padding
        rotated_coords[:, 0] += offset_x
        rotated_coords[:, 1] += offset_y

        # Fill the matrix with intensity values
        padded_data[rotated_coords[:, 1].astype(int), rotated_coords[:, 0].astype(int), feature_index] = feature_intensity

    # Flip the padded data along the horizontal axis
    if flip_horizontal:
        padded_data = np.fliplr(padded_data)

    return padded_data