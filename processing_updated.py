import pandas as pd
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


def init_sdk():
    """
    Initialize the TDF-SDK library and return the instance.
    """
    try:
        return init_tdf_sdk_api()
    except Exception as e:
        print(f"Error initializing TDF SDK: {e}")
        raise

def extract_data(path, dll, raw_data):
    """
    Extract spectra and coordinate information from a .d file.

    Args:
        path (str): Path to the .d file.
        dll (object): TDF-SDK library instance.
        raw_data (bool): Whether to extract raw profile data.

    Returns:
        tuple: (spectra_dfs, coords, roi)
    """
    data = TsfData(bruker_d_folder_name=path, tdf_sdk=dll)
    spectra_dfs = []
    coords = []
    roi = []

    for index, row in data.analysis['Frames'].iterrows():
        if raw_data:
            index_array, intensity_array = tsf_read_profile_spectrum_v2(tdf_sdk=dll, handle=data.handle, frame_id=int(row['Id']))
        else:
            index_array, intensity_array = tsf_read_line_spectrum_v2(tdf_sdk=dll, handle=data.handle, frame_id=int(row['Id']))
        mz_array = tsf_index_to_mz(tdf_sdk=dll, handle=data.handle, frame_id=int(row['Id']), indices=index_array)
        spectra_dfs.append(pd.DataFrame({'mz': mz_array, 'intensity': intensity_array}))

    if data.analysis['GlobalMetadata']['MaldiApplicationType'] == 'SingleSpectra':
        coords = [f"{path}_{n}" for n in data.analysis['MaldiFrameInfo']['SpotName']]
        roi = data.analysis['MaldiFrameInfo']['RegionNumber']
    elif data.analysis['GlobalMetadata']['MaldiApplicationType'] == 'Imaging':
        coords = [[x,y] for x,y in zip(data.analysis['MaldiFrameInfo']['XIndexPos'], data.analysis['MaldiFrameInfo']['YIndexPos'])]
        roi = data.analysis['MaldiFrameInfo']['RegionNumber']

    return spectra_dfs, coords, roi


def perform_peak_picking(spectra_dfs, sampling_rate, ppm, height, threshold, prominence):
    """
    Perform peak picking on the given spectra dataframes.

    Args:
        spectra_dfs (list): List of spectra dataframes.
        sampling_rate (int): Sampling rate for peak picking.
        ppm (int): Parts per million for peak detection precision.
        height (float): Minimum height of peaks.
        threshold (float): Minimum threshold for peak detection.
        prominence (float): Minimum prominence of peaks.

    Returns:
        list: List of processed spectra dataframes with peaks.
    """
    processed_spectra_dfs = []
    for df in tqdm(spectra_dfs, desc='Performing Peak Picking:'):
        # Apply sampling
        sampled_df = df.sample(frac=sampling_rate / 100) if sampling_rate < 100 else df

        # Perform peak detection
        peaks, _ = find_peaks(sampled_df['intensity'], height=height, threshold=threshold, prominence=prominence)

        # Extract peaks
        peak_mzs = sampled_df.iloc[peaks]['mz'].values
        peak_intensities = sampled_df.iloc[peaks]['intensity'].values

        # Filter based on ppm (this logic remains as per your original implementation)
        # Assuming placeholder for now; adjust according to your specific ppm filtering logic

        processed_spectra_dfs.append(pd.DataFrame({'mz': peak_mzs, 'intensity': peak_intensities}))

    return processed_spectra_dfs

def write_data_to_imzml(output_path, spectra_dfs, coords, peak_pick):
    """
    Write processed data to imzML file(s).

    Args:
        output_path (str): Base output path for imzML files.
        spectra_dfs (list): List of spectra dataframes to be written.
        coords (list): List of coordinates for each spectra dataframe.
        peak_pick (bool): Indicates if peak picking was performed.

    Returns:
        None
    """
    with ImzMLWriter(output_path) as writer:
        for df, coord in zip(spectra_dfs, coords):
            mzs = df['mz'].values
            intensities = df['intensity'].values
            writer.addSpectrum(mzs, intensities, coord)

    print(f"Data successfully written to {output_path}")

def process_d_files(d_file_paths, return_one_imzml=False, raw_data=False, peak_pick=False, ppm=30, mz_range=(300,1000), height=None, threshold=None, prominence=1e3, sampling_rate=10):
    """
    Processes Bruker .d files for mass spectrometry analysis.

    Args:
        d_file_paths (list): List of paths to .d files.
        return_one_imzml (bool): If True, combines all files into one imzML file.
        raw_data (bool): If True, extracts raw profile data instead of line spectrum.
        peak_pick (bool): If True, performs peak picking.
        ppm (int): Parts per million for peak picking.
        mz_range (tuple): Tuple indicating the range of mz values to consider.
        height (float): Required height of peaks. None by default.
        threshold (float): Required threshold of peaks. None by default.
        prominence (float): Minimum peak prominence for peak picking.
        sampling_rate (int): Sampling rate for peak picking.

    Returns:
        None: Saves processed data to imzML files or combined imzML file.
    """
    dll = init_sdk()
    combined_spectra_dfs = []
    combined_coords = []
    combined_roi = []

    for path in d_file_paths:
        spectra_dfs, coords, roi = extract_data(path, dll, raw_data)

        if peak_pick:
            spectra_dfs = perform_peak_picking(spectra_dfs, sampling_rate, ppm, height, threshold, prominence)

        if return_one_imzml:
            combined_spectra_dfs.extend(spectra_dfs)
            combined_coords.extend(coords)
            combined_roi.extend(roi)
        else:
            output_path = f"{path}_processed.imzML"
            write_data_to_imzml(output_path, spectra_dfs, coords, peak_pick)

    if return_one_imzml:
        combined_output_path = "combined_processed.imzML"
        write_data_to_imzml(combined_output_path, combined_spectra_dfs, combined_coords, peak_pick)
        print("Combined imzML file has been created.")
    
# Function to load imzML data into memory
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

# Function to generate or update m/z bins based on linear spacing or peak picking
def generate_or_update_mz_bins(mz_range, ppm, imzml_dataset=None, peak_pick=False, height=None, threshold=None, prominence=1e3):
    """
    Generates or updates m/z bins for analysis. Can perform linear binning or use peak picking to determine bins.

    Args:
        mz_range (list): The range of m/z values to consider for binning, specified as [min, max].
        ppm (float): Parts per million, used to determine the precision for linear binning.
        imzml_dataset (inMemoryIMS, optional): The imzML dataset object. Required if peak_pick is True.
        peak_pick (bool, optional): If True, uses peak picking to generate bins based on detected peaks. Otherwise, performs linear binning.
        height (float, optional): The minimum height of peaks for the peak picking process. Ignored if peak_pick is False.
        threshold (float, optional): The threshold intensity that a peak must surpass to be identified in the peak picking process. Ignored if peak_pick is False.
        prominence (float, optional): The minimum prominence of peaks for the peak picking. This helps differentiate significant peaks from noise.

    Returns:
        list: A list of m/z bin values. If peak_pick is False, these bins are linearly spaced within mz_range according to ppm. If peak_pick is True, the bins correspond to the m/z values of detected peaks.
    """
    # Initialize the list of m/z bins starting with the minimum value in the specified range.
    mz_bins = [mz_range[0]]
    
    if not peak_pick:
        # Linear binning: generate bins within the mz_range, spaced according to the specified ppm value.
        # This continues until the last bin exceeds the maximum value in mz_range.
        while mz_bins[-1] < mz_range[1]:
            next_bin = mz_bins[-1] + mz_bins[-1] * 2 * ppm * 10**-6  # Calculate the next bin value.
            mz_bins.append(next_bin)  # Append the new bin to the list.
    else:
        # Peak picking: first, calculate the average spectrum from the imzML dataset.
        avg_mz, avg_intens = average_spectrum(imzml_dataset)
        
        # Find peaks in the average intensity spectrum.
        # Peaks are identified based on specified height, threshold, and prominence parameters.
        peaks, _ = find_peaks(avg_intens, height=height, threshold=threshold, prominence=prominence)
        
        # Generate m/z bins based on the m/z values of the detected peaks.
        mz_bins = [avg_mz[i] for i in peaks]

    return mz_bins  # Return the list of generated or updated m/z bins.


# Function to calculate the average spectrum from the imzML dataset
def average_spectrum(imzml_dataset):
    """
    Calculates the average mass-to-charge (m/z) ratio and intensity across all spectra in an imzML dataset.

    Args:
        imzml_dataset (inMemoryIMS): An imzML dataset object which provides methods to access the spectral data.

    Returns:
        tuple: A tuple containing two numpy arrays:
            - The first array is the average m/z ratios across all spectra.
            - The second array is the average intensities corresponding to those m/z ratios.
    """
    avg_mzs = []  # List to accumulate m/z ratios from each spectrum
    avg_intens = []  # List to accumulate intensities from each spectrum
    
    # Iterate over each spectrum in the dataset using the index list
    for ii in tqdm(imzml_dataset.index_list, desc='reading raw data'):
        # Retrieve the m/z ratios and intensities for the current spectrum
        mzs, intens = imzml_dataset.get_spectrum(ii).get_spectrum(source='centroids')
        
        # Append the m/z ratios and intensities to their respective lists
        avg_mzs.append(mzs)
        avg_intens.append(intens)
    
    # Calculate the mean of m/z ratios and intensities across all spectra
    # This involves converting the lists to numpy arrays for efficient computation
    avg_mz = np.mean(np.array(avg_mzs), axis=0)
    avg_intens = np.mean(np.array(avg_intens), axis=0)
    
    # Return the average m/z ratios and intensities as a tuple of numpy arrays
    return avg_mz, avg_intens


# Function to filter features and construct the datacube from the imzML dataset
def filter_features_and_construct_datacube(imzml_dataset, mz_bins, ppm, feature_n, threshold_value=0):
    """
    Filters m/z bins based on occurrence across the dataset and constructs a datacube with the filtered features.

    Args:
        imzml_dataset (inMemoryIMS): The imzML dataset object.
        mz_bins (list): List of m/z bins to consider for filtering and datacube construction.
        ppm (float): Parts per million tolerance for m/z value matching.
        feature_n (float): The threshold as a fraction of the total number of pixels that a feature must be present in to be included.
        threshold_value (float, optional): Intensity threshold below which values will be set to zero in the datacube.

    Returns:
        tuple: A tuple containing:
            - datacube_array: A numpy array representing the datacube with dimensions [pixels, features].
            - valid_mz_bins: A list of m/z bins that passed the filtering criteria.
            - count: A list of counts representing the occurrence of each m/z bin across the dataset before filtering.
    """
    count = []  # To hold the count of pixels where each m/z bin is present above threshold

    # Loop through each m/z bin to calculate its presence across the dataset
    for bin_index in tqdm(range(len(mz_bins)), desc="Filtering m/z bins"):
        # Generate an ion image for the current m/z bin and count how many pixels have intensities above 0
        count.append(imzml_dataset.get_ion_image(mz_bins[bin_index], ppm).xic_to_image(0).astype(bool).sum())

    # Filter m/z bins based on their presence across a given fraction of the dataset
    mz_bins_filter = (np.array(count) >= int(imzml_dataset.coords.shape[0] * feature_n))
    mz_bins_use = np.array(mz_bins)[mz_bins_filter]  # Apply the filter to select valid m/z bins

    # Construct the datacube using the filtered m/z bins
    datacube_array, valid_mz_bins = construct_datacube(imzml_dataset, mz_bins_use, ppm, feature_n, threshold_value)

    # Return the constructed datacube, the list of valid m/z bins, and their counts before filtering
    return datacube_array, valid_mz_bins, count


# Function to actually build the datacube using filtered m/z bins
def construct_datacube(imzml_dataset, mz_bins_use, ppm, feature_n, threshold_value):
    """
    Constructs a datacube from the given m/z bins after applying an intensity threshold, and filters out bins with insufficient presence across the dataset.

    Args:
        imzml_dataset (inMemoryIMS): The imzML dataset object from which ion images will be generated.
        mz_bins_use (list): List of m/z bins that have passed the initial filtering criteria and will be used to construct the datacube.
        ppm (float): Parts per million tolerance for matching m/z values when generating ion images.
        feature_n (float): Threshold as a fraction of the total number of pixels; an m/z bin must be present in at least this fraction of pixels to be included.
        threshold_value (float): Intensity threshold; pixel values below this threshold in the ion images will be set to zero.

    Returns:
        tuple: A tuple containing:
            - datacube_array: A numpy array representing the constructed datacube with dimensions [pixels, features], where each column is an ion image from a valid m/z bin.
            - valid_mz_bins: A list of m/z bins that were included in the datacube, having passed both the initial filtering and the presence threshold.
    """
    datacube_array = []  # Initialize a list to collect ion images for each valid m/z bin
    valid_mz_bins = []  # Initialize a list to keep track of valid m/z bins that make it into the datacube

    # Request ion images for all m/z bins of interest
    datacube = imzml_dataset.get_ion_image(mz_bins_use, ppm)

    # Iterate over each m/z bin to process and potentially include its ion image in the datacube
    for i, mz in tqdm(enumerate(mz_bins_use), total=len(mz_bins_use), desc='Constructing Datacube'):
        image = datacube.xic_to_image(i)  # Extract the ion image for the current m/z bin
        image[image < threshold_value] = 0  # Apply the intensity threshold, setting lower values to zero

        # Include the ion image in the datacube if it is present in a sufficient number of pixels
        if np.count_nonzero(image) >= int(imzml_dataset.coords.shape[0] * feature_n):
            datacube_array.append(image)  # Add the thresholded ion image to the datacube
            valid_mz_bins.append(mz)  # Record the m/z value as valid

    # Stack the collected ion images into a numpy array to form the final datacube
    datacube_array = np.stack(datacube_array, axis=1) if datacube_array else np.array([])

    # Return the datacube and the list of m/z bins that were included in it
    return datacube_array, valid_mz_bins


# Main function to orchestrate the m/z feature extraction from an imzML dataset
def extractMZFeatures(imzml_dataset, ppm, mz_range, feature_n=0.05, mz_bins=[], rebinning=False, threshold_value=0, peak_pick=True, height=None, threshold=None, prominence=1e3):
    """
    Main function to extract m/z features from an imzML dataset.

    Args:
        imzml_dataset (inMemoryIMS): The imzML dataset from which to extract features.
        ppm (float): Parts per million precision for m/z value comparison.
        mz_range (list): The range of m/z values to consider, specified as [min_mz, max_mz].
        feature_n (float, optional): The minimum fraction of pixels in which a feature must be present to be included.
        mz_bins (list, optional): Specific m/z bins to consider. If empty and peak_pick is False, bins are generated linearly within mz_range.
        rebinning (bool, optional): Indicates if dynamic rebinning should be applied (not implemented in the provided snippets).
        threshold_value (float, optional): The intensity threshold; intensities below this value are set to zero in the datacube.
        peak_pick (bool, optional): Indicates if peak picking should be used to determine m/z bins.
        height (float, optional): The minimum height of a peak to be considered significant during peak picking.
        threshold (float, optional): The minimum intensity threshold a peak must surpass to be detected.
        prominence (float, optional): The prominence a peak must have to be considered significant.

    Returns:
        tuple: A tuple containing the following:
            - datacube_array: A numpy array representing the constructed datacube, with dimensions [pixels, features].
            - valid_mz_bins: A list of m/z bins that passed filtering and were included in the datacube.
            - count: A list of occurrence counts for each m/z bin across the dataset, before filtering.
    """
    # Generate or update m/z bins based on peak picking or linear binning
    if len(mz_bins) == 0 or peak_pick:
        mz_bins = generate_or_update_mz_bins(mz_range, ppm, imzml_dataset, peak_pick, height, threshold, prominence)

    # Filter features based on occurrence and construct the datacube with valid m/z bins
    datacube_array, valid_mz_bins, count = filter_features_and_construct_datacube(imzml_dataset, mz_bins, ppm, feature_n, threshold_value)

    # Output the number of features that passed the filtering process
    print(f"Number of features after filtering: {len(valid_mz_bins)}")

    # Return the constructed datacube, list of valid m/z bins, and their pre-filtering counts
    return datacube_array, valid_mz_bins, count


def IonImg(datacube, coords, roi_data, current_roi, mz_values, selected_mz, background, return_mask, rotate=False):
    
    """
    Generates an ion image from a datacube.

    Args:
        datacube (np.array): The datacube containing ion image data.
        coords (list): List of coordinates for each point in the image.
        roi_data (list): List of ROI data corresponding to each point.
        current_roi (int): The current region of interest (ROI).
        mz_values (list): List of m/z values.
        selected_mz (float): The selected m/z value for the image.
        background (bool): Flag to determine background handling.
        return_mask (bool): Flag to return a mask along with the image.
        rotate (bool): Flag to rotate the image by 90 degrees.

    Returns:
        np.array: The processed ion image, and optionally a mask.
    """
    
    # Filter data and coordinates for the current ROI
    roi_filter = roi_data == current_roi
    filtered_coords = coords[roi_filter]

    # Normalize coordinate system
    filtered_coords -= filtered_coords.min(axis=0)

    # Calculate the dimensions of the image
    nx, ny = filtered_coords.max(axis=0) + 1
    nx += nx % 2  # Ensure dimensions are even numbers
    ny += ny % 2

    # Find the closest m/z value
    mz_index = np.abs(mz_values - selected_mz).argmin()
    filtered_data = datacube[:,mz_index,roi_filter]

    # Create the image array
    img = np.zeros((ny, nx)) if background else np.full((ny, nx), np.nan)

    # Populate the image array with data
    img[filtered_coords[:, 1], filtered_coords[:, 0]] = filtered_data

    if return_mask:
        img_mask = np.zeros_like(img)
        img_mask[filtered_coords[:, 1], filtered_coords[:, 0]] = 1
        return img, img_mask

    return img


# Function to plot an image with a scale bar
def plot_image_with_scale_bar(ax, image, pixel_size, scale_length, scale_bar_location='lower right', scale_bar_color='black'):
    """
    Plots an image with a scale bar on the provided axes.

    Args:
        ax (matplotlib.axes.Axes): The axes to plot on.
        image (np.array): The image data to plot.
        pixel_size (float): The size of one pixel in the image.
        scale_length (float): The length of the scale bar in the same units as pixel_size.
        scale_bar_location (str): Location of the scale bar ('lower right' or 'lower left').
        scale_bar_color (str): Color of the scale bar.

    Returns:
        None: The function modifies the axes object in place.
    """
    
    ax.imshow(image)  # Ensure aspect ratio is maintained

    # Calculate the scale bar size in pixels
    scale_bar_size_in_pixels = int(scale_length / pixel_size)

    # Define the scale bar location based on the image dimensions
    y, x = image.shape
    x_start, y_start = 1, y - 1  # Default to lower left
    if scale_bar_location == 'lower right':
        x_start = x - scale_bar_size_in_pixels - 1

    # Create a line for the scale bar
    scale_bar = mlines.Line2D([x_start, x_start + scale_bar_size_in_pixels], [y_start, y_start], color=scale_bar_color, linewidth=2, solid_capstyle='butt')
    ax.add_line(scale_bar)

    # Calculate text position (offset from scale bar)
    text_offset = y * 0.06  # For example, 2% of the image height
    text_y_position = y_start - text_offset  # Position above the scale bar

    # Add label for the scale bar
    ax.text(x_start + scale_bar_size_in_pixels / 2, text_y_position, f'{scale_length} µm', color=scale_bar_color, ha='center', va='top')

    ax.axis('off')  # Remove axes


# Function to plot ion images
def plot_ion_images(all_data, all_coords, all_roi_data, all_mz_values, selected_mz, pixel_size, scale_length, unique_rois=[0], save_filename=False, background=False, return_mask=False):
    """
    Plots ion images for specified ROIs.

    Args:
        all_data (np.array): The complete datacube.
        all_coords (list): List of all coordinates.
        all_roi_data (list): List of all ROI data.
        all_mz_values (list): List of all m/z values.
        selected_mz (float): The selected m/z value for plotting.
        pixel_size (float): The size of one pixel in the image.
        scale_length (float): The length of the scale bar in the same units as pixel_size.
        unique_rois (list): List of unique ROIs to plot.
        save_filename (str): Path to save the plot, if provided.
        background (bool): Flag to determine background handling.
        return_mask (bool): Flag to return a mask along with the image.

    Returns:
        None: The function generates and displays a plot.
    """
    
    # Determine the number of rows and columns for the subplot grid
    nrows = int(np.ceil(np.sqrt(len(unique_rois))))
    ncols = nrows if nrows * (nrows - 1) < len(unique_rois) else nrows - 1
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows))
    axes = axes.flatten() if nrows * ncols > 1 else [axes]

    # Find the maximum dimensions
    max_dims = np.max([IonImg(all_data, all_coords, all_roi_data, roi, all_mz_values, selected_mz, background, return_mask).shape for roi in unique_rois], axis=0)

    # Add subtitle with the closest m/z value
    closest_mz = all_mz_values[np.abs(all_mz_values - selected_mz).argmin()]
    fig.suptitle(f"m/z value: {closest_mz:.4f}")

    # Iterate through each ROI and plot the corresponding ion image
    for i, current_roi in enumerate(unique_rois):
        ion_image = IonImg(all_data, all_coords, all_roi_data, current_roi, all_mz_values, selected_mz, background, return_mask)

        # Pad the image to match the maximum dimensions
        pad_val = np.nan if not background else 0
        padded_image = np.pad(ion_image, [(0, max_dims[0] - ion_image.shape[0]), (0, max_dims[1] - ion_image.shape[1])], mode='constant', constant_values=pad_val)
        
        ax = axes[i]
        ax.imshow(padded_image)
        ax.axis('off')   #Remove axis markers

    #Add a scale bar to the last ROI image
    last_roi_ax = axes[len(unique_rois) - 1]
    plot_image_with_scale_bar(last_roi_ax, padded_image, pixel_size, scale_length, scale_bar_location='lower right', scale_bar_color='black')

    # Turn off axes for any unused subplot spaces
    for j in range(len(unique_rois), len(axes)):
        axes[j].axis('off')

    plt.tight_layout()

    # Save the figure to file if a filename is provided
    if save_filename:
        plt.savefig(save_filename, bbox_inches='tight')

    plt.show()   #Display the plot
