# M2C2

Multiscale MALDI-2 with Contrastive Computing

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

These dependencies are automatically installed when setting up the Conda environment as described above.

## Usage

After setting up the Conda environment and activating it, you can run the scripts provided in the repository to process and analyze your mass spectrometry imaging data.
