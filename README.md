# Multimodal MALDI-2
Multimodal MALDI-2 Imaging for Spatial Neurolipidomics in Alzheimer's disease

## Installation

This project uses [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) to manage its environment and dependencies. Please ensure you have Conda installed before proceeding with the setup.

### Setting Up the Conda Environment

1. **Clone the Repository**

    Clone this repository to your local machine with Git (or download):

    ```bash
    git clone https://github.com/timtk1/M2C2.git
    cd M2C2
    ```

2. **Create the Conda Environment**

    A `environment.yml` file is provided in the root of the project directory. Use this file to create a new Conda environment with all the required dependencies:

    ```bash
    conda env create -f environment.yml
    ```

    This command reads the `environment.yml` file in the current directory and creates a new environment named as specified in the file (e.g., `myenv`).

3. **Activate the Conda Environment**

    Once the environment is created, you can activate it using:

    ```bash
    conda activate myenv
    ```

    Make sure to replace `myenv` with the actual name of your environment if it's different.

## Dependencies

This project makes use of the following key dependencies:

- **pyTDFSDK**: A Python package for processing timsTOF mass spectrometry data. It is installed directly from its GitHub repository.
- **pyimzML**: A library for reading and writing imzML files in Python, useful for mass spectrometry imaging data.
- **pyImagingMSpec**: An imaging mass spectrometry analysis package that provides tools for analyzing mass spectrometry imaging data. Note that this package was developed in Python 2, so be sure to use the package as provided in this repo, with necessary compatability changes.

## Data analysis scripts

| Script | Description |
| --- | --- |
| [Register_Fluorescence_to_MSI.m](https://github.com/timtk1/MALDI2Alz/blob/main/Register_Fluorescence_to_MSI.m) |Register optical images to MSI data, extract lipids from plaque regions |
| [RegionalAnalysis_ManualAnnotation.m](https://github.com/timtk1/MALDI2Alz/blob/main/RegionalAnalysis_ManualAnnotation.m) | Extract lipid profiles from specific brain regions |
| [PCA_MSI.m ](https://github.com/timtk1/MALDI2Alz/blob/main/Functions/%20PCA_MSI.m) | Function to fun PCA on MSI data. Used in registration script to generate MSI image with sufficient morphology for registration. |
| [ProcessMSIData.py](https://github.com/timtk1/MALDI2Alz/blob/main/RegionalAnalysis_ManualAnnotation.m) | Parse raw Bruker .d data, peak pick, bin, and export as .mat file |
| [Extract_ExpressionProfiles_BMC.py](https://github.com/timtk1/MALDI2Alz/blob/main/RegionalAnalysis_ManualAnnotation.m) |Extract expression profiles from [Lee et al., BMC Genomics.](https://link.springer.com/article/10.1186/s12864-024-10434-8)  |






 ## Example data
Data can be accessed [here](https://databank.illinois.edu/datasets/IDB-4907703).
- `Animal_1_5xFAD_s1.mat`: A MATLAB array of 50 micron spatial resolution imaging of whole brain slice from a 5xFAD animal.
- `mz_bins_use_neg.mat`: A MATLAB array of the m/z channels all MSI images (whole brain slice, 50 micron spatial resolution) were binned to in order to enable comparison
- `Animal3_S18_HR.mat`: A MATLAB array of high-spatial-resolution (5 micron) imaging of a 5xFAD mouse hippocampus and cortex. Due to the large dataset, 22 m/z channels ('mz_features_22.mat') are included.
- `Animal5_S18_HR.mat`: A MATLAB array of high-spatial-resolution (5 micron) imaging of a WT mouse hippocampus and cortex. Due to the large dataset, 22 m/z channels ('mz_features_22.mat') are included.
- `mz_features_22.mat`: A MATLAB array of the 22 m/z channels included in the high spatial resolution imaging data


